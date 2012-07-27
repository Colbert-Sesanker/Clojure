;; Colbert Sesanker 2012 
;; a scheme interpreter in clojure in response to Peter Norvig's Lispy: http://norvig.com/lispy.html

(use 'clojure.string)

(defn append [val vec]
  (assoc vec (count vec) val))

(defn cdr [vec]
  (subvec vec 1))
 
(def globals (ref {"+" +, "-" -, "*" *,"/" /, "not" not,
      ">" >, "<" <, ">=" >=, "<=" <=,"=" =,  "equal?" =, 
      "eq?" =, "length"count, "cons" cons, "car" first, "cdr" cdr, 
      "append" append, "list" list, "list?" vector?, "null?" nil?, "symbol?" symbol?}))

(defn add [var-val obj]
  (dosync (ref-set obj (apply assoc @obj var-val)))
    obj)

(defn zip-vec [key-s val-s] 
  (if (vector? key-s)       
    (vec (interleave key-s val-s))
    [key-s val-s]))

(defn eval 
  ([x] (eval x globals))
  ([x env]
   (cond 
     (number? x)  x
     (not (vector? x)) (@env x)    ;;looks up x in env map                   
     (= (x 0) "qoute") (next x)
     (= (x 0) "if")    (eval (if (eval (x 1), env) (x 2) (x 3)), env)
     (= (x 0) "set!")  (add [(x 1) (eval (x 2) env)] env)
     (= (x 0) "define")(add [(x 1) (eval (x 2) env)] env)
     (= (x 0) "lambda")(fn [& args] (eval (x 2), (add (zip-vec (x 1) args) env)))
     (= (x 0) "begin") (last (map #(eval % env) (next x)))
     :else  (let [exps (vec (map #(eval % env) x)) fun (exps 0)]
              (apply fun (next exps))))))

(defn tokenize [s]
   (split 
     (trim 
       (replace 
         (replace s #"\)" " ) " ) 
       #"\(" " ( "))
   #"\s+"))        

(defn len [v] ;; calculates length of parsed s-expression  
    (let [cnt #(if (vector? %) (+ 2 (count %)) 1) ]
      (inc (reduce +  
            (map cnt v)))))    
 
;; this is a recursive mess, clean-up if you can
(defn rcat [tokens l fun]  ;;l is a reference to the parsed output vector     
  (if (not= (tokens 0) ")")  
    (if (= (tokens 0) "(") 
      (let [a (fun tokens) tok-len (inc (len a))]
          (rcat (subvec tokens tok-len), (append a l),  fun))   
        (rcat (cdr tokens), (append (fun tokens) l),  fun)) 
     l))  
                                                 
(defn read-from [tokens]
  (assert (not= 0 (count tokens)))
    (let [token (tokens 0)]  
     (cond  
       (= "(" token)  (rcat (cdr tokens) [] read-from)                 
       ( = ")" token) (throw (Exception. "unexpected ( "))
       :else (atomic token))))
     
(defn atomic [token]   
   (try (Integer/parseInt token) 
    (catch Exception e 
      (try (Float/parseFloat token) 
        (catch Exception e
           token)))))  

(defn to-string [exp]
  (if (vector? exp)
    (str "(" (join " " (map to-string exp)) ")")
      (str exp))) 

(defn scheme-repl []
  (print (format "%s: " "min-lisp->> "))
  (flush)
  (let [a (read-line) parse #(read-from (tokenize %))] 
  (println (to-string (eval (parse a))))    
     (if (not= a "exit")       
        (recur)
        (println "bye"))))    
                      
;(scheme-repl)
;(define area (lambda (r) (* 3.141592653 (* r r))))
;(area 3)