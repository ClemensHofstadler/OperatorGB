# coding: utf-8
"""
NCPoly
================

A class that provides basic funcionality for noncommutative polynomials

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

from .NCMonomial import NCMonomial

###############################################################################
# NCPoly
###############################################################################
class NCPoly:
    def __init__(self,f,*args):
        if args:
            self.lt = f
            self.tail = args[0]
            if len(args) > 1:
                self.cofactors = args[1]
        else:
            d = f.monomial_coefficients()
            if d:
                monomials = [NCMonomial(d[key],key) for key in d.keys()]
                monomials = sorted(monomials)
            else:
                monomials = [NCMonomial(0,'*',intern=True)]
            self.lt = monomials[-1]
            self.tail = monomials[:-1]
            
            self.cofactors = None
        self.redundant = False
############################################################################
    def copy(self):
        return NCPoly(self.lt.copy(),[mon.copy() for mon in self.tail])
############################################################################
    def zero():
        return NCPoly(NCMonomial.zero(),[],[])
############################################################################
    def __repr__(self):
        s = "(" + str(self.lt) + ", ["
        for m in self.tail:
            s += str(m) + ", "
        if self.tail:
            s = s[:-2]
        s += "])"
        return s
############################################################################
    def __eq__(self,other):
        if other is None:
            return False
        if len(self.tail) != len(other.tail):
            return False
        return self.lt == other.lt and all([m1 == m2 for (m1,m2) in zip(self.tail,other.tail)])
############################################################################
    def __hash__(self):
        return hash((self.lt,) + tuple(self.tail))
############################################################################
    def makeMonic(self,i = None):
        """
        Really update self
        """
        lc = self.lt.coeff
        if lc != 1:
            self.lt.coeff = 1
            for m in self.tail:
                m.coeff /= lc
        if i is not None:
            self.cofactors = [(NCMonomial(1/lc,"*",intern=True),i,NCMonomial.one())]
############################################################################
    def lmul(self,m):
        """
        Only update copy
        """
        f = self.copy()
        f.lt.mon = m[:-1] + f.lt.mon
        for mon in f.tail:
            mon.mon = m[:-1] + mon.mon
   
        return f
############################################################################
    def rmul(self,m):
        """
        Only update copy
        """
        f = self.copy()
        f.lt.mon += m[1:]
        for mon in f.tail:
            mon.mon += m[1:]
   
        return f
############################################################################
    def lrmul(self,l,r):
        f = self.copy()
        f.lt.mon = l[:-1] + f.lt.mon + r[1:]
        for mon in f.tail:
            mon.mon = l[:-1] + mon.mon + r[1:]
        
        return f
############################################################################
    def toNormal(self,Parent=None):
        return self.lt.toNormal(Parent) + sum([_.toNormal(Parent) for _ in self.tail])
