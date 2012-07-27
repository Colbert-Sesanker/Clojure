#Example Output: 
    user=> (scheme-repl)
    min-lisp->> : (define area (lambda (r) (* 3.141592653 (* r r))))
    clojure.lang.Ref@992fa5
    min-lisp->> : (area 3)
    28.274334
    min-lisp->> : (define fact (lambda (n) (if (<= n 1) 1 (* n (fact (- n 1))))))
    clojure.lang.Ref@992fa5
    min-lisp->> : (fact 10)
    8800
    min-lisp->> : exit
    bye
    nil
    user=> 


