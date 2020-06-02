# OperatorGB
Mathematica and Sage package to enumerate non-commutative Groebner bases together with tracing of cofactors.
To get started in Mathematica see the file 'Documentation.nb' and to use the Sage version of the package look at the 'Documentation.py' file. 

Version 1.2 includes (for both versions of the package):
  * Different monomial orderings (DegLex, MultiLex and in Mathematica also weighted DegLex)
  * An interreduction procedure for polynomials
  * Groebner basis computations together with tracing of cofactors
  * Chain criterion to delete redudant ambiguities
  * Quivers (not necessarily with unique labels)
  * Computing the signature of polynomials with respect to a quiver 
  * Checking compatibility as well uniform compatibliity of polynomials with a quiver
  * A 'Certify' command that does everything from checking compatibility, to computing a Groebner basis as well as a cofactor representation all fully automatically.
  * Q-completion procedure (so far only in Mathematica) (new in this version)
-----------
<b>NEW in Version 1.2.1:</b>

Mathematica:
 * IntegerCoeffQ function
 * CheckCertificate function

<b>NEW in Version 1.2.2:</b>

Mathematica:
 * Info option has been changed to global variable VerboseOperatorGB (now 3 possibilities - no info, some info, all info)
 * Intersection of (two-sided) ideal with another (two-sided) ideal or subalgebra
 * Can automatically find elements in an ideal to which *-cancellability is applicable
 * Reorganised stuff
 
----------- 
TODO:

Both:
 * Option to compute reduced GB 
 * Test file

Mathematica:
  * Trace cofactors for the intersections
  
Sage:
  * Q-order-compatibility
  * Q-completion
  * IntegerCoeffQ function
  * Info -> VerboseOperatorGB
  * Intersections
  * *-cancellability
