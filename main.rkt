 #lang racket
(require "declarations.rkt")
(require "drawing-routine.rkt")
(require "testcases.rkt")

(define (calcF pd cd d mp mc)
  (if (= d 0) 0 (/ (* g mp mc (- cd pd)) (* d d d))))

(define (calcPoint tree p)
  (let*([cx (vec-x (gnode-posn tree))]
        [cy (vec-y (gnode-posn tree))]
        [px (vec-x (particle-posn p))]
        [py (vec-y (particle-posn p))]
        [mp (particle-mass p)]
        [mc (gnode-mass tree)]
        [d (sqrt (+ (* (- cx px) (- cx px)) (* (- cy py) (- cy py))))]
        [rx (calcF px cx d mp mc)]
        [ry (calcF py cy d mp mc)])
    (vec rx ry)))

(define (calcInd side tree p)
  (cond ((null? tree) (vec 0 0))
      ((null? (gnode-subtrees tree)) (calcPoint tree p))
      (else (let*([cx (vec-x (gnode-posn tree))]
        [cy (vec-y (gnode-posn tree))]
        [px (vec-x (particle-posn p))]
        [py (vec-y (particle-posn p))]
        [mp (particle-mass p)]
        [mc (gnode-mass tree)]
        [d (sqrt (+ (* (- cx px) (- cx px)) (* (- cy py) (- cy py))))]
        [rat (/ d side)])
    (cond ((> rat theta)
           (let*([fx (calcF px cx d mp mc)]
                 [fy (calcF py cy d mp mc)])
             (vec fx fy)))
          (else (let*([rx (sum (lc (vec-x (calcInd (/ side 2) sub p)) : sub <- (gnode-subtrees tree)))]
                      [ry (sum (lc (vec-y (calcInd (/ side 2) sub p)) : sub <- (gnode-subtrees tree)))])
                  (vec rx ry))))))))

(define (range? p area)
  (if (and (and (< (vec-x (particle-posn p)) (bbox-rux area)) (>= (vec-x (particle-posn p)) (bbox-llx area)))
           (and (< (vec-y (particle-posn p)) (bbox-ruy area)) (>= (vec-y (particle-posn p)) (bbox-lly area))))
      #t #f))

(define (buildTree area particles)
  (let*([l (lc x : x <- particles @(range? x area))]
        [lbx (bbox-llx area)]
        [ubx (bbox-rux area)]
        [lby (bbox-lly area)]
        [uby (bbox-ruy area)]
        [mx (/ (+ lbx ubx) 2)]
        [my (/ (+ lby uby) 2)]
        [tmass (sum (lc (particle-mass x) : x <- l))]
        [cx (if (> tmass 0) (/ (sum (lc (* (particle-mass p) (vec-x (particle-posn p))) : p <- l )) tmass) mx)]
        [cy (if (> tmass 0) (/ (sum (lc (* (particle-mass p) (vec-y (particle-posn p))) : p <- l )) tmass) my)])
    (cond ((null? l) '())
          ((singleton l) (gnode (particle-mass (car l)) (particle-posn (car l)) '()))
          (else
           (let*([t1 (buildTree (bbox lbx lby mx my) l)]
                 [t2 (buildTree (bbox mx lby ubx my) l)]
                 [t3 (buildTree (bbox lbx my mx uby) l)]
                 [t4 (buildTree (bbox mx my ubx uby) l)]
                 [tfin (lc nod : nod <- (list t3 t4 t1 t2) @(not (null? nod)))])
             (gnode tmass (vec cx cy) tfin))))))

(define (calcForces area tree particles)
  (lc (calcInd (- (bbox-rux area) (bbox-llx area)) tree p) : p <- particles))

(define (effect-force p f)
  (let*([m (particle-mass p)]
        [fx (vec-x f)]
        [fy (vec-y f)]
        [ax (/ fx m)]
        [ay (/ fy m)]
        [sx (vec-x (particle-posn p))]
        [sy (vec-y (particle-posn p))]
        [ux (vec-x (particle-velocity p))]
        [uy (vec-y (particle-velocity p))]
        [vx (+ ux (* ax timeslice))]
        [vy (+ uy (* ay timeslice))]
        [dx (+ sx (* ux timeslice) (* 0.5 ax (* timeslice timeslice)))]
        [dy (+ sy (* uy timeslice) (* 0.5 ay (* timeslice timeslice)))])
    (particle m (vec dx dy) (vec vx vy))))

(define (moveparticles particles forces)
  (zipwith effect-force particles forces))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; singlestep :: (Particle) -> (Particle)


(define (singlestep particles)
  (let* ([initialArea (bounding-box particles)]
         [tree (buildTree initialArea particles)] 
         [forces (calcForces initialArea tree particles)])
    (moveparticles particles forces)))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;main will send the current vectors for printing
;calculate new vectors and recurse

;;main :: (Particle) -> Action
(define (main ps)
  (define (main-helper ps i)
    (cond [(> i iter) (display "Done")]
          [else (let*
                    ([ps-next (singlestep ps)])
                  (if (= (remainder iter drawtime) 0)
                      (begin 
                        (draw-particles ps)
                        (main-helper ps-next (+ i 1)))
                      (main-helper ps-next (+ i 1))))]))
  (main-helper ps 0))

(main testList2)
