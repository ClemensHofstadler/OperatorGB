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
TODO:

Mathematica:
  * Info: version that prints only some information
  * Certify/Groebner without Info => many new lines. Get rid of this
  
Sage:
  * Q-order-compatibility
  * Q-completion
