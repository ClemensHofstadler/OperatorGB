# coding: utf-8
"""
F4
================

TODO Description
AUTHOR:

- Clemens Hofstadler (2017-08-08)

"""

#############################################################################
#  Copyright (C) 2019 Clemens Hofstadler (clemens.hofstadler@liwest.at).    #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################


############################################################################
# PhD Examples
############################################################################
#F = FreeAlgebra(ZZ,2,['a','b'])
#(a,b) = F.gens()
#G5 = [a**2-1,b**5-1,(a*b*a*b**2)**2 -1]
#SetUpRing(b,a)
#G = F4(G5,4)
#G6 = [a*a-1,b^5-1,(a*b*a*b*a*b^4)^2-1]

############################################################################
# Big Hartwig
############################################################################
vars = ['a','adj_a','aa','adj_aa','b','adj_b','bb','adj_bb','c','adj_c','cc','adj_cc','x','adj_x','u','adj_u','v','adj_v']
F = FreeAlgebra(ZZ,len(vars),vars)
(a,adj_a,aa,adj_aa,b,adj_b,bb,adj_bb,c,adj_c,cc,adj_cc,x,adj_x,u,adj_u,v,adj_v) = F.gens()
p = aa*a*b*c*cc
q = c*cc*bb*aa*a
m = a*b*c
adj_m = adj_c*adj_b*adj_a
adj_p = adj_cc*adj_c*adj_b*adj_a*adj_aa
adj_q = adj_a*adj_aa*adj_bb*adj_cc*adj_c
pinvA = [a*aa*a-a,aa*a*aa-aa,adj_aa*adj_a-a*aa,adj_a*adj_aa-aa*a]
pinvB = [b*bb*b-b,bb*b*bb-bb,adj_bb*adj_b-b*bb,adj_b*adj_bb-bb*b]
pinvC = [c*cc*c-c,cc*c*cc-cc,adj_cc*adj_c-c*cc,adj_c*adj_cc-cc*c]
rest = [m*x*m-m,x*m*x-x,adj_x*adj_m-m*x,adj_m*adj_x-x*m,p*q*p*q-p*q,adj_q-adj_a*a*p*v,c*adj_c*adj_p-q*u]
rest2 = [-adj_a+adj_a*adj_aa*adj_a,-adj_aa+adj_aa*adj_a*adj_aa,
        -adj_b+adj_b*adj_bb*adj_b,-adj_bb+adj_bb*adj_b*adj_bb,
        -adj_c+adj_c*adj_cc*adj_c,-adj_cc+adj_cc*adj_c*adj_cc,
        -adj_m+adj_m*adj_x*adj_m,-adj_x+adj_x*adj_m*adj_x,
        -adj_q*adj_p+adj_q*adj_p*adj_q*adj_p,q-adj_v*adj_p*adj_a*a,adj_u*adj_q+p*c*adj_c]
I = pinvA + pinvB + pinvC + rest + rest2
claims = [p*q*p-p,
        q*p*q-q,
        adj_a*a*p*q-adj_a*adj_aa*adj_bb*adj_cc*adj_c*adj_cc*adj_c*adj_b*adj_a*adj_aa*adj_a*a,
        q*p*c*adj_c - c*adj_c*adj_cc*adj_c*adj_b*adj_a*adj_aa*adj_a*adj_aa*adj_bb*adj_cc*adj_c
        ]
SetUpRing(F.gens())
#SetUpRing(F.gens()[:-4],F.gens()[-4:])
