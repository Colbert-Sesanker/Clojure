;; Author: Colbert Sesanker
;; email: sesanker0@gmail.com
;; Multi-classifier SVM with 21 classes (21 amino-acids) 
;; Test examples are peptide sequences with gaps represented by X
;; 
;; Start a REPL, load in name-space and functions. 
;; The main functions are:

;; 1. (com-pare "A" "B" "C" "D")  
;;    Enter the three arguments as paths to plain text files with peptide sequences 
;; Inputs:
;;    "A" : Path to protein used to train the classifier.  
;;    "B" : Path to protein with gaps to test.    
;;    "C" : Path to the same protein as argument "B" without any gaps. Used to score the accuracy of the classifier
;; Output:
;;    "D" : Name of output file. Outputs result in the current directory.
;;
;;     Try: (com-pare "order.txt" "orderX.txt" "order.txt")  
;;          Generate these files using the procedure below. Notice the files "A" and "C" are the same.
;;          
;;          This means the file used to train the classifier is the same as the one used to score it.
;;          Since this data is artificially generated, the model should find patterns in the data, and it does.
;;          To generate this data in the REPL type (order-protein 100 ["LRKDC" "YDEST" "ESGPI"] "order.txt")                     
;;          The order-protein function shuffles the three 5-letter chunks, "LRKDC" "YDEST" "ESGPI", concatenates them,
;;          shuffles them again, concatenates them to the previous and so on 100 times. This yields a string of
;;          1500 characters written out as plain text to the file order.txt, in the current directory.
;;          To generate orderX.txt simply replace any 4 Gs with X in order.txt using a find/replace in any text editor.
;;          
;;    To understand why the classifier works so well on the order.txt data it is necessary to understand 
;;    the classification process: When com-pare takes in a training peptide sequence, "A", it breaks it up into
;;    5 letter chunks, e.g., LRKDC, also called  neighborhoods. It then partitions all of the chunks into disjoint
;;    classes based on the amino acid in the center, yielding 21 classes. For example LRKDC belongs to the K-class. 
;;    In general, XX?XX belongs to the ?-class, where X is any amino acid. 
;;    Each ?-class is a set of training examples that trains the ?-model to guess ? from the four surrounding
;;    amino acids. For example, an A-class could be {KLASG MNART} where KLASG and MNART are examples for A.
;;    The features for each example AB?CD, where A,B,C,D do are NOT amino acid abbreviations, but are any amino acids.
;;    are: AD, BC, A, D, B, C, AB and CD. Since any of A, B, C and D can have 21 values, they are
;;    (21^2)*4 + 21*4 = 1848 possible features for each example. The hash-map, protein-neighborhood, associates
;;    each of the possible features with a unique index using the following language for the features:

;;    (1)                    AD, BC, A, D, B, C, AB, CD = AD02, BC01, 1A, D2, 1B, C1, 1AB, CD1

;;    The features of the hash-map are in the form on the RHS of (1) and two element combinations are computed 
;;    using the cartesian product in clojure.contrib.combinatorics. The hash-map, protein-neighborhood
;;    has 1848 key value pairs in the form {AL01 33, 2D 66, 1LL 1330, R2 1847, ...}, where the pairs are unordered. 
;;    The a sparse feature vector in R^1848 is constructed from each example. The operations are defined
;;    so only non zero values need to appear. A feature vector could be
;;    {66 1, 44 1, 1010 1, 33 1, 3939 1, 0 1, 99 1, 2333 1}. They are always 8 non-zero entries in each feature vector
;;    since 8 possible features can be extracted for each example. Notice the values are always 1, atypical of SVMs.

;;    This section requires a basic understanding of Support Vector Machines. The ?-model is trained by projecting
;;    onto (via the inner product) feature vectors derived from the ?-class, and away from feature vectors derived from
;;    the 20 other classes. Clearly they're a lot more negative examples than positive examples for each class. 
;;    From the training data, "A", com-pare builds a model for each of the 21 amino acids and tests each model on "B".
;;    Thus multiple amino acids may be proposed for each gap, X. com-pare returns a sequence, in file "D", with N+1 
;;    elements. N is the number of X's in "B". The first element is always the overall score, correct/total, 
;;    how many X's were guessed correctly based on file "C". The rest of the ordered elements correspond to the Xs in file "B"; 
;;    the ith element of the sequence in "D" corresponds to the ith X in file "B", indexed from left to right. Create an example
;;    output "D" by typing: (com-pare "order.txt" "orderX.txt" "order.txt" "out.txt") in the REPL. out.txt reads:

;;                                 ([:score 1] [\G "TG"] [\G "TG"] [\G "TG"] [\G "TG"])

;;    For each element after the score in the form [a b], a corresponds to the actual value of the gap and b is an 
;;    ordered sting, ordered from RIGHT TO LEFT of guesses for the gap; In the case of [\G "TG"], \G means G is the
;;    actual value and "TG" means G is the most likely guess and T is the second likely guess. 
;;                

;; 2. show-plots: graphs learning statistics for a given amino acid using Incanter
;;    try  (show-plots "A" "Anole.txt") and compare to (show-plots "A" "order.txt")
;;          after making ordered data via (order-protein 100 ["LRKDC" "YDEST" "ESGPI"] "order.txt")  
;;    
;; The vector of models, model-matrix, (a vector of vectors) transforms the binary classification problem into 
;; a 21 classifier problem. Notice the training of the models is intrinsically parallel.

;; A number of functions are provided for manipulating protein text data in the supporting functions section. 

;; Overall, This project does exactly what it's designed to, but real protein sequences are not uniformly
;; structured locally in five piece chunks. Nonetheless, it introduces a general way to define features and
;; many of the ideas and code are reusable to solve other problems.


(ns AminoAcidPredictor     
   (:use clojure.contrib.combinatorics)
   (:use [clojure.contrib.generic.functor :only (fmap)])  
   (:use [incanter.charts                 :only (scatter-plot add-lines)])  
   (:use [incanter.core                   :only (view)])
   (:import (java.io FileReader BufferedReader))
   (:import (org.apache.commons.lang StringUtils)))  			 
       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  Static Definitions  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def amino-acids       ;; 21 amino acids including "U"
  ["R" "H" "K" "D" "E" "S" "T" "N" "Q" "C" "U" "G" "P" "A" "V" "I" "L" "M" "F" "Y" "W"])  

(def protein-neighbors ;; Sequence of all two letter words from the alphabet "amino-acids" 
	 (cartesian-product amino-acids amino-acids))

;; Creates a hash-map that associates a unique integer to each possible feature
(def protein-neighborhood (let [clump-2   (fn [[a b]] (str a b))  s (map clump-2 protein-neighbors)] 
  (zipmap (reduce into (fmap vec [ (map str s (repeat "01")) (map str s (repeat "02")) 
                                   (map str (repeat "1") s) (map str s (repeat "1")) 
                                   (map str amino-acids (repeat "1")) (map str amino-acids (repeat "2"))
                                   (map str  (repeat "1") amino-acids ) (map str (repeat "2") amino-acids) ]))                           
                       (range))))   
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   Supporting functions   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; reports function evaluation time in ms
(defn bind-time [f & args]
  "function -> args -> evluation time"
  (let [start (. System (nanoTime))         
        j     (apply f args)
        end   (. System (nanoTime))]
                        (/ (- end start) 1000000.0)))

;; writes string f to a file (without quotes)
(defn out-file [f file-name]                
  {:pre [(string? f) (string? file-name)]}
  (binding [*out* (java.io.FileWriter. file-name)]                       
          (prn (symbol f))))

;; reads file and outputs contents as a string
(defn read-file [file-path]
 "File Path (string) - > contents (string)"
 {:pre [(string? name)]}
 (.trim (reduce str (line-seq (BufferedReader. 
                                (FileReader. file-path))))))

;; randomly repalaces n amino-acids in a protein string with "X" 
(defn x-aa [protein n]                       
 "protein (string) -> number of X replacements (int) -> protein with Xs (string)" 
  (reduce str (apply assoc (vec protein) 
                      (vec (interleave (repeatedly n #(rand-int (count protein))) (repeat n \X))))))

;; reads in a protein file then writes out same-file with animo-acids 
;; randomly replaced with Xs
(defn noise [protein-file proteinX-file n]   
   (let [protein (read-file protein-file)]   
     (out-file (x-aa protein n) proteinX-file))) 

;; Do not call by itself, will not terminate.
(defn rand-aa "Returns a random amino" []     
  (repeatedly #(rand-nth amino-acids)))

;; creates a random protein with "size" amino acids
;; writes to file with name file-out
(defn rand-protein [size file-out] 
 {:pre [(integer? size)]}          
  (out-file (apply str (take size (rand-aa))) file-out))

;; shuffles neighborhoods to always associate 
;; the same features about amino-acids 
(defn order-protein [repeats seed-vector file-out]   
 {:pre [(integer? repeats) (vector? seed-vector)]}   
      (out-file (apply str (repeatedly (int repeats) #(apply str (shuffle seed-vector)))) file-out))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; End supporting functions;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Below are  5 math functions from Mark Reid: http://mark.reid.name/sap/online-learning-in-clojure.html
;; Sparse vectors are represented by maps with key value pairs. The key represents the index in the
;; sparse vector, and the value is the value of the index. Note, only nonzero values are stored. 

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; End Math functions;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Below is a scheme to transform peptide sequences (FASTA format) into a training set that can 
;; train an SVM to fill gaps in peptide sequences

(defn feature-index    ;; Get the unique index for a given feature
   "String -> Int"     ;; uses hash-map protein neighborhood as a function
   [protein-neighbor] 
   (protein-neighborhood protein-neighbor))
        
 
(defn neighborhood-seq
  "Peptide Sequence (String) -> Target amino-acid -> neighborhood length -> Neighborhoods "
  [protein target radius]   ;; returns length 5 neighborhood about target          
  (let [ r                     (apply str (repeat radius '.))
         target-regex          (reduce str [r "[" target "]" r])
	       neighborhood-sequence (re-seq (re-pattern target-regex) protein)]
         neighborhood-sequence))

(defn make [[a b _ c d]]       ; extracts features from the neighborhood  
 "match (String) -> tuple of neighborhood features"
 [(str 2 a) (str 1 b) (str 1 a b) (str b c 0 1 ) (str a d 0 2) (str c d 1) (str c 1) (str d 2)])

(defn index-tuple [[a b c d e f g h]]     
"feature vector (vector) -> sparse vector of local features extracted about an amino-acid (vector)"
 {(feature-index a) 1, (feature-index b) 1, (feature-index c) 1, (feature-index d) 1, 
  (feature-index e) 1, (feature-index f) 1, (feature-index g) 1 ,(feature-index h) 1})

(defn sparse-vec-seq [neighborhood-seq]   ;; returns sequence of sparse vectors
  "neighborhood-sequence (seq string) -> sparse-vector-sequence (seq {}) "
 (map (comp index-tuple make) neighborhood-seq))

(defn target-vectors [amino-acid protein]
  "amino-acid (char) -> protein (string) -> sparse-vector-sequence (seq {})"
 (sparse-vec-seq (neighborhood-seq protein amino-acid 2)))

(defn error-vectors [amino-acid protein]
   "amino-acid (char) -> protein (string) -> sparse-vector-sequence (seq {})"
 (sparse-vec-seq (neighborhood-seq protein (str "^" amino-acid) 2)))

(defn pos-example [sparse-vector]
 "Sparse vector, a map of key value pairs (map)-> Map assigning sparse vector to positive class (Map)"
  {:y  1 :x sparse-vector})

(defn neg-example [sparse-vector]
 "Sparse vector, a map of key value pairs (map)-> Map assigning sparse vector to positive class (Map)"
  {:y -1 :x sparse-vector})
 
(defn examples [target-vectors error-vectors] ;; A sequence of examples used for training the model in {:y sgn :x sparse vector} format
 (into (map pos-example target-vectors) (map neg-example error-vectors))) ;; examples are unordered and independent of each other

;; Given a file path reads a peptide sequence text file and turns it into training data
;; or given an amino acid a file path and option returns examples or target vectors for a given amino acid
(defn protein-to-examples 
 ([protein-file-path]  
  (let [ protein        (read-file protein-file-path)
         targets-matrix (into [] (map #(target-vectors % protein) amino-acids)) ; returns matrix of examples for all amino-acids
         errors-matrix  (into [] (map #(error-vectors % protein) amino-acids))
         example-matrix (into [] (map #(examples % %2) targets-matrix errors-matrix))] ; Specify data to get out                                        
                          example-matrix)) ; vector of examples ordered in the same way as the amino acids
 
      ;; overload: if amino-acid is specified return examples for that amino-acid               
      ([amino-acid protein-file-path data]   
      {:pre [ (integer? data) (or (= 0 data) (= 1 data)) ]} 
        (let [  protein  (read-file protein-file-path)
                targets  (target-vectors amino-acid protein)          
                errors   (error-vectors amino-acid protein)      
                examples (if (= data 0) (examples targets errors)
                                         targets) ]
                          examples)))    
 
;; Projects (inner product) a model onto a feature vector and returns the amino acid corresponding to
;; the model and the projection, in a tuple, if the projection is greater than 1.
;; This projections are used to sort candidates for a given feature vector based on projection.
(defn projection-tuple  
 "amino-acid (string) -> model (map) -> {string #} or nil" 
  [amino-acid w example] 
 (let [ projection  (inner w example)]                      
	 (if (<= 1 projection) [amino-acid projection] nil)))     

;; creates a vector of amino-acid projection pairs given an amino acid
;; its corresponding model and test data
(defn ordered-solution-vector 
  "amino acid (string) -> model (map) -> test-data (string)"
  [amino-acid model test-data]
   (vec (map #(projection-tuple amino-acid model %) test-data)))  




;; Functions hinge-loss and correct are from Mark Reid. See: http://mark.reid.name/sap/online-learning-in-clojure.html 

(defn hinge-loss ;; Returns the hinge loss of the weight vector w on the given example
	"model w (map) -> Training example (map) -> non-negative number   "
	[w example] (max 0 (- 1 (* (:y example) (inner w (:x example))))))
	
(defn correct  ;; Returns a corrected version of the weight vector w
	[w example t lambda]
	(let [x   (:x example)
		  y   (:y example)
		  w1  (scale (- 1 (/ 1 t)) w)
		  eta (/ 1 (* lambda t))
		  r   (/ 1 (Math/sqrt lambda))]
		(project (add w1 (scale (* eta y) x)) r)))

(defn update
	"folding step to train the model"
	[model example]
	(let [ t        (:step   model)
	       w        (:w      model)             
	       fix      (> (hinge-loss w example) 0) ]
   
	       { :w     (if fix (correct w example t 0.0001) w),			  
	         :step     (inc t)}))

(defn train
	"Returns a model trained from the initial model on the given examples"
	[initial-model examples]
	(reduce update initial-model examples))

;; creates a vector of models (sequence of vectors) for each of 21 amino acids in 
;; the same order as amino acids defined above
(defn create-model-matrix [file-path] 
 "File-path (string) -> ordered vector of 21 models (sequence of vectors)" 
  (let [  initial-model {:step 1, :w {}} 
          example-matrix (protein-to-examples file-path)
          model-matrix   (map #(:w (train initial-model %)) example-matrix) ]          
                         model-matrix))

;; Tests the projection of each model onto the features about each X in the 
;; test data. Note, test data is a string here not a file path
(defn ord-sol-matrix [model-matrix test-data]
  "sequence of vectors -> Test data string -> vector of tuples "
  {:pre [ (= (count model-matrix) (count amino-acids)) ]}
   (vec (map #(ordered-solution-vector %  %2 test-data) amino-acids model-matrix)))         

;; returns the column vectors of an m*n matrix: 
;; ([a11 a12 ... a1n] [a21 a22 ... a2n] ...) --->  ([a11 a21 ... an1] [a12 a22 ... an2]...)
(defn mat-cols [vec-seq]       
  (let [n   (count vec-seq)   
        m   (count (first vec-seq))
        f  #(into [] (map nth (vec %) (repeat n %2)))]    
                     (vec (map #(f vec-seq %) (range m)))))

;; Sorts proposed solutions by projection and concatenates into a string
(defn structure [mat-cols]
  (map (comp #(apply str %) #(map first %) #(sort-by second %)) mat-cols))
 
;; Returns a sequence of strings. Each string is ordered by projection from right to left 
(defn ord-sol-seq [training-file-path testing-file-path]
  (let [  model-matrix     (create-model-matrix training-file-path)
          test-data        (protein-to-examples "X"  testing-file-path 1)] ;; 1 returns features from test data          
          (structure (mat-cols (ord-sol-matrix model-matrix test-data)))))  

;;x is the protein with Xs, s the original, indicies is an index bin
(defn x-seq-index [x s indicies i]        
  (if (= x ((vec s) i)) 
          (into indicies [i]) 
                   indicies))

;; returns the indicies in actual protein corresponding to the X's
(defn findX [x s]      
  (reduce #(x-seq-index x s % %2) [] (range (count s)))) 

;; Finds the ith amino acid in test data without Xs corresponding to the ith X the test data with Xs
(defn actual [t tX]
  (apply str (vec (map #((vec t) %) (findX \X tX)))))

;; Retuns the ratio of how many gaps preducted correctly to the total number of gaps
(defn score [train test]                
  (let [ p        (ord-sol-seq train test)
         t        (read-file train)
         tX       (read-file test)        
         actual   (actual t tX)          
         correct  (map #(if (= % (last %2)) 1 0) actual p)]                
                          (assert (-> correct (= 0) (not))) 
                          (/ (reduce + correct) (count correct))))
;; See introduction
(defn com-pare [t rX r out]
 (out-file (str 
             (cons [:score (score t rX)] (map #(vec [% %2]) 
                                            (actual (read-file t) (read-file rX)) (ord-sol-seq r rX)))) out))
 
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Stats;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn build-stats
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
			  :errors   (into  errors  [(if error (inc (last errors)) (last errors))]) } ))

(defn stats-model
	"Returns statistics about the trained model and the model itself"
	[initial-model examples]
	(reduce build-stats initial-model examples))

(defn list-scatter-plot ;; Makes a scatter plot from two number lists (Incanter)
  "x-list-> y-list -> x-label (String) -> y-label (String)-> ScatterPopup"
  [x-list y-list  x-label y-label]  
   (view
    (scatter-plot x-list y-list
                  :x-label x-label 
                  :y-label y-label
                  :legend true)))	
(defn stats
	"Returns a scatter plot of correction step against the percent correct "
	[model]
	(let [	      last-step     (:step model) ;;fit model here and measure curvature
		      steps         (take (:step model) (iterate inc 1))
		      features      (:features model)
		      errors        (:errors model) 
                      error-to-step (map #(- 1 (/ (float %) %2)) errors steps) ]                      
			(list-scatter-plot  steps  error-to-step "Steps" "Percentage Correct")
                        (list-scatter-plot  steps  features "Steps" "Features")
                        (list-scatter-plot  steps  errors "Steps" "Errors")))

(defn show-plots [amino-acid train-file-path]
   (let [initial-model {:lambda 0.0001, :step 1, :w {}, :features [0], :errors [0]} 
     examples (protein-to-examples amino-acid train-file-path 0)     ;; 0 returns examples to train
     model (stats-model initial-model examples)]
      (stats model)))







		




