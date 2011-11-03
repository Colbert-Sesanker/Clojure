(ns machine
  (:use clojure.contrib.combinatorics)) 


(clojure.core/refer 'clojure.core)
(ns machine
  (:use (java.io FileReader BufferedReader))         
  (:use (clojure.contrib.combinatorics) )     
  (:use (cern.colt.matrix.tfloat.impl.SparseFloatMatrix1D))
  (:use (cern.jet.math.tfloat FloatFunctions))
  (:use (edu.emory.mathcs.utils ConcurrencyUtils)))


(defn n-tuple-list [n list] "Integer- > TupleList -> List" 
 (map #(nth % n) list))

(defn stringseq [a b] "List -> List -> StringList"
 (map str a  b))

(defn stringseq-tuple [tuple-list] "TupleList -> StringList"
 (map str (map first tuple-list) (map second tuple-list)))

(def protein-neighbors 
 (let [ proteins (repeat 2 ["R" "H" "K" "D" "E" "S" "T" "N" "Q" "C" "U" "G" "P" "A" "V" "I" "L" "M" "F" "Y" "W"])]
	(cartesian-product (first proteins) (second proteins))))

(def protein-neighborhood (let [ s (stringseq-tuple protein-neighbors)] 
  (zipmap (into (stringseq s (repeat "1")) (stringseq s (repeat "2"))) 
          (iterate inc 1))))

(defn get-protein-neighbor-index [protein-neighbor] "String -> Int"
	(protein-neighborhood protein-neighbor))

(defn target-neighbors [string]
 (let [ matches (re-seq #"..[A].." string) ]
    matches))

(defn error-neighbors [string]
 (let [ matches (re-seq #"..[S].." string) ]
    matches))

(def l 
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

(def examples
 (interleave (map pos-example target-vectors) (map neg-example error-vectors)))





