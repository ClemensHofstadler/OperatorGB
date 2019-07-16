# OperatorGB
Mathematica package to enumerate non-commutative Groebner bases together with tracing of cofactors.

Changes compared to version 1.0.0:
  * Adapted information message for SetUpRing, Weight and SortedQ
  * Certify returns linear combination in terms of generators of the ideal (and not in terms of monic versions of the generators anymore)
  * Added multigraded lexicographic order (as in NCAlgebra)
  * Certify computes the Groebner basis iteratively one iteration at a time
  * Certify also provides multigraded lexicographic order (first block: all variables appearing in the claims - assumed to be the known variables; second block: the rest)
  * Bug fix in adj
  * Added interreduction procedure and integrated it in Certify (but does not really speed things up)
  ---------
  * DeleteRedundant: Delete dedundant ambiguities (Mora's approach)
  * New record for BigHartwig: Interreduced system; multilex order; 8 iterations, MaxDeg: 25, Criterion: True => 8.5 sec
