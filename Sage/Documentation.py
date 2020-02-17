# coding: utf-8
"""
Documentation
================
In this file, we prove a simple operator identity to illustrate the syntax of the OperatorGB package. The example is taken from the
Mathematica file 'Documentation.nb', which describes the Mathematica version of this package. Here, we focus solely on how to execute this
example in Sage using this package. For details on the example itself please see the 'Documentation.nb' file.

Before loading this file into Sage, the OperatorGB package has to be loaded. To this end, place the pckage in a directory
where Sage can find it (e.g. the current working directory) and type

    from OperatorGB import *
    
Then, the whole functionality provided by the package should be available in Sage and this file
can be executed. This can be done by changing to the directory where this file is located and executing

    load("Introductory_Example.py")
    
AUTHOR:

- Clemens Hofstadler (2020-02-16)

"""

#############################################################################
#  Copyright (C) 2020 Clemens Hofstadler (clemens.hofstadler@jku.at).       #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

"""
In the following, we discuss the example from the Mathematica documentation.
Note that we write 'aa' and 'bb' instead of 'a^{-}' and 'b^{-}', resepctively.

In order to be able to enter our polynomials for the assumptions and the claim, we first have to define a free algebra. Then, we can
enter the polynomials as elements of this free algebra. We use the FreeAlgebra class provided by Sage.
"""
F = FreeAlgebra(QQ,4,['a','aa','b','bb'])
(a,aa,b,bb) = F.gens()
assumptions = [a*aa*a - a, b*bb*b - b, aa*a*b*bb*aa*a*b*bb - aa*a*b*bb]
claim = a*b*bb*aa*a*b - a*b
print("\nWe have just entered the assumptions and the claims as elements of the free algebra.")
print("Our assumptions are: " + str(assumptions))
print("Our claim is %s \n" %str(claim))

"""
Next, we have to define the quiver encoding our operators. This is done by providing triples of the form
(l,s,t) where 'l' is the label, 's' is the source and 't' is the target of an edge.
"""
Q = Quiver([('a','Vv','Vw'),('aa','Vw','Vv'),('b','Vu','Vv'),('bb','Vv','Vu')])
Q.plot()

"""
We can compute the signature of polynomials and check compatibility as well as uniform compatibility of
polynoimals with our quiver Q as follows. The methods are named just as in Mathematica, however, note that
here the methods are class methods, i.e. they have to be called with the '.' operator.
"""
print("The signature of the claim is " + str(Q.QSignature(claim)))
print("\nWe can also check whether the assumptions are uniformly compatible with Q and whether the claim is compatible with Q:")
print("Uniform compatibility of the assumptions: " + str([Q.uniformlyCompatibleQ(f) for f in assumptions]))
print("Compatibility of the claim: %s \n" % str(Q.compatibleQ(claim)))

"""
Now we have to set up the monomial order for our Groebner basis computations. So far, we provide a degree lexicographic ordering,
 which is used whenever only one list of variables is passed as an argument to the 'SetUpRing' method. If several lists are passed as
 arguments to the 'SetUpRing' method, then a block ordering is chosen. For more details on these orderings please see the
 'Documentation.nb' file.
"""
print("If we define a degree lexicographic ordering, we get the following information:")
SetUpRing([a,aa,b,bb])
print("If we define a block ordering, we the this information:")
SetUpRing([a,aa],[b,bb])
print("However, for now we switch back to the degree lexicographic ordering.")
SetUpRing([a,aa,b,bb])

"""
Now we can compute a Greobner basis. To this end, we first have to define an empty list where the cofactors will be safed in.
The other arguments are the same as in Mathematica (also all optional arguments are equal).
"""
print("\nNow we can compute a Groebner basis...\n")
cofactorsG = []
G = Groebner(cofactorsG,assumptions)
print("The Groebner basis is a list of polynomials in the free algebra")
print(G)

"""
When we have a Groebner basis, we can redcue our claim using the ReducedForm method
"""
print("\nThe reduced form of the claim is ")
cofactorsF = []
ReducedForm(cofactorsF,G,claim)
print("and we obtain a 'linear combination' of the form " + str(cofactorsF))

"""
Multiplying out such a linear combination or rewriting it follows the same syntax as in Mathematica.
"""
print("Multiplying out the linear combination yields %s \n" % str(MultiplyOut(cofactorsF)))
certificate = Rewrite(cofactorsF,cofactorsG)
print("After rewriting, the linear combination looks as follows:")
print(certificate)
print("Multiplying out the rewritten linear combination still yields " + str(MultiplyOut(certificate)))

"""
The package also provides a 'Certify' method that does all the steps described above at once.
"""
print("\nWe can do all these steps at once using the Certify command...\n")
certificate = Certify(assumptions,claim,Q)
print("\nThe output of Certify is:\n" + str(certificate))

"""
We can easily check the correctness of our certificate using the method 'CheckCertificate'
"""
print("\nWe can check whether the certificate is indeed correct...\n")
CheckCertificate(certificate,claim,assumptions)
