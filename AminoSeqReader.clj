(ns AminoSeqReader
   (:use clojure.contrib.combinatorics)
   (:use [incanter.charts :only ( bar-chart scatter-plot add-lines)])  
   (:use [incanter.core :only (view)])
   (:import (java.io FileReader BufferedReader)))  					 


;; Supporting functions
(defn n-tuple-list [n list] "Integer- > TupleList -> List" 
 (map #(nth % n) list))

(defn stringseq [a b] "List -> List -> StringList"
 (map str a  b))

(defn stringseq-tuple [tuple-list] "TupleList -> StringList"
 (map str (map first tuple-list) (map second tuple-list)))


(defn make-scatter-plot ;; Makes a scatter plot from two lists of data
  "x-list-> y-list -> String -> String -> ScatterPopup"
  [x-list y-list  x-label y-label]  
   (view
    (scatter-plot x-list y-list
                  :x-label x-label 
                  :y-label y-label
                  :legend true)))	


;; Mathematical functions from Mark Reid: http://mark.reid.name/sap/online-learning-in-clojure.html
;; Will replace with incanter libraries
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

;;

;;
;;Devise scheme to transform FASTA peptide sequences into a training set that can train an SVM to fill gaps
;;in peptide sequences

(def protein-neighbors 
 (let [ proteins (repeat 2 ["R" "H" "K" "D" "E" "S" "T" "N" "Q" "C" "U" "G" "P" "A" "V" "I" "L" "M" "F" "Y" "W"])]
	(cartesian-product (first proteins) (second proteins))))

(def protein-neighborhood (let [ s (stringseq-tuple protein-neighbors)] 
  (zipmap (into (stringseq s (repeat "1")) (stringseq s (repeat "2"))) 
          (iterate inc 1))))

(defn get-protein-neighbor-index [protein-neighbor] "String -> Int"
	(protein-neighborhood protein-neighbor))

(defn target-neighbors [string] ;; returns neighborhoods about desired amino acid
 (let [ matches (re-seq #"..[A].." string) ]
    matches))

(defn error-neighbors [string] ;;returns neighborhoods about undesired amino acid
 (let [ matches (re-seq #"..[S].." string) ]
    matches))

(def l ;;sequence FASTA format for protein
"LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY")

(defn make [matches]
 [(str (nth matches 1) (nth matches 3) 1 ) (str (first matches) (last matches)  2)])

(defn index-tuple [made]
 {(get-protein-neighbor-index (first made)) 1, (get-protein-neighbor-index (second made)) 1})

(def target-features 
 (map make (target-neighbors l)))

(def error-features 
 (map make (error-neighbors l)))

(def target-vectors
 (map index-tuple target-features))

(def error-vectors
 (map index-tuple error-features))

(defn pos-example [sparse]
{:y 1 :x sparse})

(defn neg-example [sparse]
{:y -1 :x sparse})

(def examples ;; examples used for training
 (interleave (map pos-example target-vectors) (map neg-example error-vectors)))

;; Legacy Pegasos algorythym described by Mark Reid http://mark.reid.name/sap/online-learning-in-clojure.html

(defn hinge-loss
	"Returns the hinge loss of the weight vector w on the given example"
	[w example] (max 0 (- 1 (* (:y example) (inner w (:x example))))))
	
(defn correct
	"Returns a corrected version of the weight vector w"
	[w example t lambda]
	(let [x   (:x example)
		  y   (:y example)
		  w1  (scale (- 1 (/ 1 t)) w)
		  eta (/ 1 (* lambda t))
		  r   (/ 1 (Math/sqrt lambda))]
		(project (add w1 (scale (* eta y) x)) r)))

(defn update
	"Returns an updated model by taking the last model, the next training 
	 and applying the Pegasos update step"
	[model example]
	(let [lambda (:lambda model)
		  t      (:step   model)
		  w      (:w      model)
		  errors (:errors model)
		  error  (> (hinge-loss w example) 0)]
		(do 
			
			{ :w      (if error (correct w example t lambda) w), 
			  :lambda lambda, 
			  :step   (inc t), 
			  :errors (if error (inc errors) errors)} )))

(defn train
	"Returns a model trained from the initial model on the given examples"
	[initial examples]
	(reduce update initial examples))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;End Pegasos;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn stats
	"Returns a scatter plot of correction step against the correction "
	[model]
	(let [step (:step model)
		      features   (count (keys (:w model)))
		      errors (:errors model) 
                      error-to-step (- 1 (/ (float errors) step))]                      
			(make-scatter-plot  step  error-to-step "Steps" "Percentage Correct")
                        (make-scatter-plot  step  features "Steps" "Features")))

(defn show-plots []
   (let [start {:lambda 0.0001, :step 1, :w {}, :errors 0} 
     model (train start examples)]
      (stats model)))


		

