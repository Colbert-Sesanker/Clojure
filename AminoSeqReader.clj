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






