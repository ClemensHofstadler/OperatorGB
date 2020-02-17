# coding: utf-8
"""
NCMonomial
================

A class that provides basic funcionality for noncommutative monomials, i.e. words in the free monoid.

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

from sage.all import FreeAlgebra,copy,QQ,ZZ

from . import GlobalData
from .MonomialOrder import SortedQ

############################################################################
# NCMonomial
############################################################################
class NCMonomial:
    def __init__(self,coeff,monomial,intern = False):
        if str(monomial) == '1':
            self.coeff = coeff
            self.mon = '*'
        else:
            if intern:
                self.mon = monomial
                self.coeff = coeff
            else:
                mon = '*'
                for v in str(monomial).split('*'):
                    j = v.find('^')
                    if j != -1:
                        mon += (v[:j] + '*')*ZZ(v[j+1:])
                    else:
                        mon += v + '*'
  
                self.mon =  mon
                self.coeff = coeff
############################################################################
    def copy(self):
        return NCMonomial(self.coeff,copy(self.mon),intern=True)
############################################################################
    def __lt__(self,other):
        return SortedQ(self.mon,other.mon)
############################################################################
    def __repr__(self):
        return "(" + str(self.coeff) + ", " + self.mon + ")"
############################################################################
    def __eq__(self,other):
        return self.mon == other.mon and self.coeff == other.coeff
############################################################################
    def __hash__(self):
        return hash(self.mon)
############################################################################
    def __mul__(self,other):
        if other in QQ:
            self.coeff *= other
            return self
        else:
            return NCMonomial(self.coeff*other.coeff,self.mon[:-1] + other.mon,intern = True)
############################################################################
    def one():
        return NCMonomial(1,'*',intern=True)
############################################################################
    def zero():
        return NCMonomial(0,'1')
############################################################################
    def toNormal(self,Parent=None):
        s = str(self.coeff)
        s += self.mon[:-1]
        if not Parent:
            Parent = FreeAlgebra(QQ,len(GlobalData.WordOrder),GlobalData.WordOrder)
        return Parent(s)
