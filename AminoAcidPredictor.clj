;;Author: Colbert Sesanker
;;email sesanker0@gmail.com
;;Any Legacy material is explicitly commented as such


(ns AminoAcidPredictor
   (:use clojure.contrib.combinatorics)
   (:use [incanter.charts :only (scatter-plot add-lines)])  
   (:use [incanter.core :only (view)])
   (:import (java.io FileReader BufferedReader)))  			 
       
   
;; Supporting functions

(defn n-tuple-list [n list] "Integer- > TupleList -> List" 
 (map #(nth % n) list))

(defn stringseq [a b] "List -> List -> StringList"
 (map str a  b))

(defn stringseq-tuple [tuple-list] "TupleList -> StringList"
 (stringseq (map first tuple-list) (map second tuple-list)))


(defn list-scatter-plot ;; Makes a scatter plot from two number lists
  "x-list-> y-list -> x-label (String) -> y-label (String)-> ScatterPopup"
  [x-list y-list  x-label y-label]  
   (view
    (scatter-plot x-list y-list
                  :x-label x-label 
                  :y-label y-label
                  :legend true)))	

(defn bind-time [f & args]
  (let [start (. System (nanoTime)) ;;reports function time in ms
        j     (apply f args)
        end   (. System (nanoTime))]
                        (/ (- end start) 1000000.0)))


;; The following are Mathematical functions from Mark Reid: http://mark.reid.name/sap/online-learning-in-clojure.html
;; These will be replaced with  incanter matrix libraries for increased speed
;; Sparse vectors are represented are maps with keys value pairs coresponding index value pairs in the vector

(defn add
	"Returns the sparse sum of two sparse vectors x y"
	[x y] (merge-with + x y))

(defn inner
	"Computes the inner product of the sparse vectors (hashes) x and y"
	[x y] (reduce + (map #(* (get x % 0) (get y % 0)) (keys y))))

(defn norm
	"Returns the l_2 norm of the (sparse) vector v"
	[v] (Math/sqrt (inner v v)))

(defn scale
	"Returns the scalar product of the sparse vector v by the scalar a"
	[a v] (zipmap (keys v) (map * (vals v) (repeat a))))

(defn project
  "Returns the projection of a parameter vector w onto the ball of radius r"
	[w r] (scale (min (/ r (norm w)) 1) w))

;; End of functions from Mark Reid

;; Devise scheme to transform FASTA peptide sequences into a training set that can 
;; train an SVM to fill gaps in peptide sequences

(def amino-acids 
  ["R" "H" "K" "D" "E" "S" "T" "N" "Q" "C" "U" "G" "P" "A" "V" "I" "L" "M" "F" "Y" "W"])  

(comment
(defn make-random-control-data [it alphabet]  ;; recursive function to make random string of length it*cardinality(alphabet) from a list of the alphabet's letters
  {:pre [(integer? it) (not (neg? it))]} ;;returns stack overflow for large it use loop recur
  (if (= it 1) 
      (reduce str (shuffle alphabet))
      (str (reduce str (shuffle alphabet)) 
            (make-random-control-data (dec it) alphabet))))
)

(defn make-random-control-data [num alphabet]
  (apply str (repeatedly num #(apply str (shuffle alphabet)))))

(def random-control-data
     (make-random-control-data (int 1000) amino-acids))

(defn make-ordered-control-data [i seed]
    (reduce str (into [] (repeat i seed))))

(def ordered-control-data 
    (make-ordered-control-data (int 10000) "KIALK" ))  

(def protein-neighbors ;;  
	 (cartesian-product amino-acids amino-acids))

(def protein-neighborhood (let [ s (stringseq-tuple protein-neighbors)] 
  (zipmap (into (stringseq s (repeat "1")) (stringseq s (repeat "2"))) 
          (iterate inc 1))))

(defn get-protein-neighbor-index [protein-neighbor] "String -> Int"
  (protein-neighborhood protein-neighbor))

(defn target-neighbors 
  "Peptide Sequce (String) -> Neighborhoods "
  [string] ;; returns length 5 neighborhood about Alanine
  (let [ target-matches (re-seq #"..[A].." string) ]
    target-matches))

(defn error-neighbors 
  "Peptide Sequce (String) -> Neighborhoods "
  [string]  ;;returns length 5 neighborhoods about #{Amino-Acids}\{A}
  (let [ error-matches (re-seq #"..[^A].." string) ]
    error-matches))

;;
(comment (def l ;;sequence FASTA format for protein
"
"))

(def l ordered-control-data )

(defn make [matches] 
 "sequence of neighborhoods (Sring-seq) -> length 2 sequence of features in neighboorhod"
 [(str (nth matches 1) (nth matches 3) 1 ) (str (first matches) (last matches)  2)])

(defn index-tuple [made]
 "sequence of neighborhoods (Sring-seq) -> length 2 sequence of features in neighboorhod"
 {(get-protein-neighbor-index (first made)) 1, (get-protein-neighbor-index (second made)) 1})

(def target-features 
;;Creates a unordered sequence of features for each target neighborhood
 (map make (target-neighbors l)))

(def error-features 
 (map make (error-neighbors l)))

(def target-vectors
 (map index-tuple target-features))

(def error-vectors
 (map index-tuple error-features))

(defn pos-example [sparse-vector]
 "Sparse vector, a map of key value pairs (map)-> Map assigning sparse vector to positive class (Map)"
  {:y  1 :x sparse-vector})

(defn neg-example [sparse-vector]
 "Sparse vector, a map of key value pairs (map)-> Map assigning sparse vector to positive class (Map)"
  {:y -1 :x sparse-vector})

(def examples ;; A map of examples used for training the model in {:y sgn :x sparse vector} format
 (into (map pos-example target-vectors) (map neg-example error-vectors)))

;; The following is SVM algorythym based on one described by Mark Reid 
;;http://mark.reid.name/sap/online-learning-in-clojure.html 
;;Functions hinge-loss, correct and train are from Mark's implementation

(defn hinge-loss ;; Returns the hinge loss of the weight vector w on the given example
	"model w (map) -> Training example (map) -> non-negative number   "
	[w example] (max 0 (- 1 (* (:y example) (inner w (:x example))))))
	
(defn correct  ;;Returns a corrected version of the weight vector w
	[w example t lambda]
	(let [x   (:x example)
		  y   (:y example)
		  w1  (scale (- 1 (/ 1 t)) w)
		  eta (/ 1 (* lambda t))
		  r   (/ 1 (Math/sqrt lambda))]
		(project (add w1 (scale (* eta y) x)) r)))

(defn update
	"Returns an updated model by taking the last model, the next training 
	 and applying the Pegasos update step. Allocates errors and features indexed by each time step"
	[model example]
	(let [      lambda   (:lambda model)
		    t        (:step   model)
		    w        (:w      model)
                    features (count (keys (:w model)))
		    errors   (:errors model)
		    error    (> (hinge-loss w example) 0)]
   
			{ :w        (if error (correct w example t lambda) w), 
			  :lambda   lambda, 
			  :step     (inc t), 
                          :features (into (:features model) [features])
			  :errors   (into  errors  [(if error (inc ( last errors)) (last errors))]) } ))

(defn train
	"Returns a model trained from the initial model on the given examples"
	[initial-model examples]
	(reduce update initial-model examples))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn stats
	"Returns a scatter plot of correction step against the percent correct "
	[model]
	(let [	last-step     (:step model)
		      steps         (take (:step model) (iterate inc 1))
		      features      (:features model)
		      errors        (:errors model) 
                      error-to-step (map #(- 1 (/ (float %) %2)) errors steps) ]                      
			(list-scatter-plot  steps  error-to-step "Steps" "Percentage Correct")
                        (list-scatter-plot  steps  features "Steps" "Features")
                        (list-scatter-plot  steps  errors "Steps" "Errors")))

(defn show-plots []
   (let [start {:lambda 0.0001, :step 1, :w {}, :features [0], :errors [0]} 
     model (train start examples)]
      (stats model)))


		


