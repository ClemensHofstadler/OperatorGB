# OperatorGB
Mathematica package to enumerate non-commutative Groebner bases together with tracing of cofactors.

Changes in version 1.1.0 compared to version 1.0.0:
  * Adapted information message for SetUpRing, Weight and SortedQ
  * Certify returns linear combination in terms of generators of the ideal (and not in terms of monic versions of the generators anymore)
  * Added multigraded lexicographic order (as in NCAlgebra)
  * Certify computes the Groebner basis iteratively one iteration at a time
  * Certify also provides multigraded lexicographic order (first block: all variables appearing in the claims - assumed to be the known variables; second block: the rest)
  * Bug fix in adj
  * Added interreduction procedure and integrated it in Certify.
  * DeleteRedundant: Delete dedundant ambiguities (Mora's approach)
  * Got rid of complicated data structure ReductionSystem and cleaned code up
  * Just 4 fun: implemented (terribly slow) F4 that can trace cofactors
  ---------
Changes compared to version 1.1.0:
  * Changed order in Overlap ambiguitiy from (ABC,C,A,gi,gj) to (ABC,A,C,gi,gj)
  * Certify now only rewrites the cofactors actually appearing in the certificate
  --------
Changes since last meeting:
  * Mathematica:
 * Sage:
    - Implemented Quiver
    - Implemented Faugere's algorithm to compute RRef
  ________
TODO:
  * Compute minimal Groebner basis (Thm. 5.3.10 in PhD thesis)
  * Discuss F5 criteria
  * Chain criterion
  * Mathematica:
    - adapt Documentation
  
  
