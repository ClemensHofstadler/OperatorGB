# coding: utf-8
"""
NCMonomial
================

A class that provides basic funcionality for noncommutative monomials, i.e. words in the free monoid.

AUTHOR:

- Clemens Hofstadler (2019-07-08)

"""

#############################################################################
#  Copyright (C) 2019 Clemens Hofstadler (clemens.hofstadler@liwest.at).    #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import absolute_import
from sage.all import FreeAlgebra,copy,QQ

import GlobalData
from MonomialOrder import SortedQ

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
                        mon += str(v) + '*'
  
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
    def toNormal(self):
        s = str(self.coeff)
        s += self.mon[:-1]
        F = FreeAlgebra(QQ,len(GlobalData.WordOrder),GlobalData.WordOrder)
        return F(s)
