# coding: utf-8
"""
Introductory example
================
In this file we prove a simple operator identity to illustrate the syntax of the OperatorGB package. The example is taken from the
Mathematica file 'Documentation.nb', which describes the Mathematica version of this package. Here, we focus solely on how to execute this
example in Sage using this package. For details on the example itself please see the 'Documentation.nb' file.

Before loading this file into Sage the OperatorGB package has to be loaded. To this end, place the 'src' folder in a directory
where Sage can find it (or simply change to the directory 'src' using the 'cd' command) and type

    from OperatorGB import *
    
Then, the whole functionality provided by the package should be available in Sage and this file
can be executed. This can be done by changing to the directory where this file is located
(again by using the 'cd' command) and executing

    load("./Introductory_Example.py")
    
AUTHOR:

- Clemens Hofstadler (2019-10-17)

"""

#############################################################################
#  Copyright (C) 2019 Clemens Hofstadler (clemens.hofstadler@liwest.at).    #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

"""
The syntax for this package is basically the same as the syntax for the equally named Mathematica
package. However, there are some small differences:
    - in Mathematica we have
            'g = ReducedForm[vars,G,exp]'
        whereas here we have
            '(g,vars) = ReducedForm(G,exp)'
    - in Mathematica we have have
            '(signatureAssumptions,signatureClaims,normalForms,certificate) = Certify(assumptions,claims,Q)'
        whereas here we have
            'certificate = Certify(assumptions,claims,Q)' if all claims can be reduced to zero. Otherwise 'False' is returned.
    - in Mathematica the first argument of Groebner just has to be a name. In a list with this name the cofactors will then be safed. Here, we explicitly have to pass a list.


In the following we discuss the example from the Mathematica documentation.
Note that we write 'aa' and 'bb' instead of 'a^{-}' and 'b^{-}', resepctively.

As a first step we have to define the quiver encoding our operators. This is done by providing triples of the form
(l,s,t) where 'l' is the label, 's' is the source and 't' is the target of an edge.
"""
Q = Quiver((('a','V','W'),('aa','W','V'),('b','U','V'),('bb','V','U')))
print "We have just defined the following quiver Q:\n" + str(Q)
print "\nWe can also plot the quiver..."
Q.plot()

"""
In order to be able to enter our polynomials for the assumptions and the claim, we first have to define a free algebra. Then we can
enter the polynomials as elements of this free algebra. Here, we use the FreeAlgebra class provided by Sage.
"""
F = FreeAlgebra(QQ,len(Q.vars()),Q.vars())
(a,aa,b,bb) = F.gens()
assumptions = [a*aa*a - a, b*bb*b - b, aa*a*b*bb*aa*a*b*bb - aa*a*b*bb]
claim = a*b*bb*aa*a*b - a*b
print "\nWe have just entered the assumptions and the claims as elements of the free algebra."
print "Our assumptions are: " + str(assumptions)
print "Our claim is %s \n" %str(claim)

"""
We can compute the signature of polynomials and check compatibility as well as uniform compatibility of
polynoimals with our quiver Q as follows. The methods are named just as in Mathematica, however, note that
here the methods are class methods, i.e. they have to be called with the '.' operator.
"""
print "The signature of the claim is " + str(Q.QSignature(claim))
print "\nWe can also check whether the assumptions are uniformly compatible with Q and whether the claim is compatible with Q:"
print "Uniform compatibility of the assumptions: " + str([Q.uniformlyCompatibleQ(f) for f in assumptions])
print "Compatibility of the claim: %s \n" % str(Q.compatibleQ(claim))

"""
Now we have to set up the monomial order for our Groebner basis computations. So far, we provide a degree lexicographic ordering,
 which is used whenever only one list of variables is passed as an argument to the 'SetUpRing' method. If two lists are passed as
 arguments to the 'SetUpRing' method, then a block ordering is chosen. For more details on these orderings please see the
 'Documentation.nb' file.
"""
print "If we define a degree lexicographic ordering, we get the following information:"
SetUpRing(Q.vars())
print "If we define a block ordering, we the this information:"
SetUpRing(Q.vars()[:2],Q.vars()[2:])
print "However, for now we switch back to the degree lexicographic ordering."
SetUpRing(Q.vars())

"""
Now we can compute a Greobner basis. To this end, we first have to define an empty list where the cofactors will be safed in.
The other arguments are the same as in Mathematica (also all optional arguments are equal.
"""
print "\nNow we can compute a Groebner basis...\n"
cofactors = []
G = Groebner(cofactors,assumptions)
print "The Groebner basis is a list of polynomials in the free algebra"
print G

"""
When we have a Groebner basis, we can redcue our claim. As already described above, the syntax of the method 'ReducedForm' is
slightly different as in Mathematica.
"""
print "\nUsing this Groebner basis we can reduce to claim..."
(g,vars) = ReducedForm(G,claim)
print "\nThe reduced form of the claim is " + str(g)
print "and we obtain a 'linear combination' of the form " + str(vars)

"""
Multiplying out such a linear combination or rewriting it follows the same syntax as in Mathematica.
"""
print "Multiplying out the linear combination yields %s \n" % str(MultiplyOut(vars))
linearComb = Rewrite(vars,cofactors)
print "After rewriting, the linear combination looks as follows:"
print(linearComb)
print "Multiplying out the rewritten linear combination still yields " + str(MultiplyOut(linearComb))

"""
As the Mathematica package, also this package provides a 'Certify' method that does all the steps described above at once. However, as noted
above the output of this method is different than in Mathematica.
"""
print "\nWe can do all these steps at once using the Certify command...\n"
certificate = Certify(assumptions,claim,Q)
print "\nThe output of Certify is:\n" + str(certificate)

"""
We can easily check the correctness of our certificate using the method 'CheckCertificate'
"""
print "\nWe can check whether the certificate is indeed correct...\n"
CheckCertificate(certificate,claim,assumptions)
