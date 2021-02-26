# OperatorGB
Mathematica and Sage package to enumerate non-commutative Groebner bases together with tracing of cofactors.
To get started in Mathematica see the file 'Documentation.nb' and to use the Sage version of the package look at the 'Documentation.py' file. 

# Mathematica version of the package:

Version 1.3 includes:
  * Different monomial orderings (DegLex, MultiLex and in Mathematica also weighted DegLex)
  * An interreduction procedure for polynomials
  * Groebner basis computations together with tracing of cofactors
  * Chain criterion to delete redudant ambiguities
  * Quivers (not necessarily with unique labels)
  * Computing the signature of polynomials with respect to a quiver 
  * Checking compatibility as well uniform compatibliity of polynomials with a quiver
  * A 'Certify' command that does everything from checking compatibility, to computing a Groebner basis as well as a cofactor representation all fully automatically.
  * Q-completion procedure
  * Several procedures to do ideal exploration (Elimination orderings, ideal intersections, subalgebra intersection)
  * (Bi-)module Groebner basis computations
  * Simple-to-use interface for applying cancellability
  * Simple-to-use interface for diagram chases. 

# Sage version of the package:

Version 1.2 includes:
  * Different monomial orderings (DegLex, MultiLex and in Mathematica also weighted DegLex)
  * An interreduction procedure for polynomials
  * Groebner basis computations together with tracing of cofactors
  * Chain criterion to delete redudant ambiguities
  * Quivers (not necessarily with unique labels)
  * Computing the signature of polynomials with respect to a quiver 
  * Checking compatibility as well uniform compatibliity of polynomials with a quiver
  * A 'Certify' command that does everything from checking compatibility, to computing a Groebner basis as well as a cofactor representation all fully automatically.
  
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
  * Module Gröbner bases
  * Ideal exploration
  * Interfaces for properties of operators/diagram chases
