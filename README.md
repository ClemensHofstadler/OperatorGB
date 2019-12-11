# OperatorGB
Mathematica and Sage packages to enumerate non-commutative Groebner bases together with tracing of cofactors.
To get started in Mathematica see the file 'Documentation.nb' and to use the Sage version of the package look at the 'Introductory_Example.py' file. 

Version 1.1.0 includes (for both programs):
  * Different monomial orderings (DegLex, MultiLex and in Mathematica also weighted DegLex)
  * An interreduction procedure for polynomials
  * Groebner basis computations together with tracing of cofactors
  * Some deletion criteria (Mora's approach + chain criterion)
  * Quivers (not necessarily with unique labels)
  * Computing the signature of polynomials with respect to a quiver 
  * Checking compatibility as well uniform compatibliity of polynomials with a quiver
  * A 'Certify' command that does everything from checking compatibility, to computing a Groebner basis as well as a cofactor representation all fully automatically.

DONE:

For both:
  * Make documentation modular
  
 Mathematica: 
  * change output for Certify depending on 'Info'
  * add functionality to autmatically generate trivial quiver from given polynomials

TODO:

Mathematica:
  * Info: version that prints only some information
  * Certify only reduce to 0 once. (Maybe try again and save shortest linear comb.)
  * Certify/Groebner without Info => many new lines. Get rid of this
  * Update doc in sourcecode.
  
Sage:
  * same as above
  * change output for Certify depending on 'Info'
  * clean up package
