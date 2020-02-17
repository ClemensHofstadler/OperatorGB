# coding: utf-8
"""
Ambiguites
================

Module that provides functionality to compute ambiguities of words

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

import itertools
from time import time

from .MonomialOrder import SortedQ,flatten

############################################################################
# ambiguities
############################################################################
class Overlap:
    """
    fiC - Afj
    """
    def __init__(self,ABC,A,C,i,j):
        self.ABC = ABC
        self.A = A
        self.C = C
        self.i = i
        self.j = j
        self.deg = len(ABC.split('*'))
        self.min = min(i,j)

    def __repr__(self):
        return "Overlap(" + str(self.ABC) + ", " + str(self.A) + ", " + str(self.C) + ", (" + str(self.i) + ", " + str(self.j) + "))"
############################################################################
class Inclusion:
    """
    fi - AfjC
    """
    def __init__(self,ABC,A,C,i,j):
        self.ABC = ABC
        self.A = A
        self.C = C
        self.i = i
        self.j = j
        self.deg = len(ABC.split('*'))
        self.min = min(i,j)

    def __repr__(self):
        return "Inclusion(" + str(self.ABC) + ", " + str(self.A) + ", " + str(self.C) + ", (" + str(self.i) + ", " + str(self.j) + "))"
############################################################################
def GenerateAmbiguities(words, newPart = None, MaxDeg = -1):
    if newPart != None:
        return __generateAmbiguities__(words,newPart,MaxDeg)
    
    amb = [generateOverlaps(v,w) for (v,w) in itertools.product(words,repeat=2)]
    amb += [list(generateInclusions(v,w)) for (v,w) in itertools.product(words,repeat=2)]

    amb = flatten(amb)
   
    if MaxDeg > 0:
        amb = [a for a in amb if a.deg-2 <= MaxDeg]

    return amb
############################################################################
def __generateAmbiguities__(oldPart,newPart,MaxDeg):
    
    amb = [generateOverlaps(v,w) + generateOverlaps(w,v) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [list(generateInclusions(v,w)) + list(generateInclusions(w,v)) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [generateOverlaps(v,w) for (v,w) in itertools.product(newPart,repeat=2)]
    amb += [list(generateInclusions(v,w)) for (v,w) in itertools.product(newPart,repeat=2)]
    
    amb = flatten(amb)
    
    if MaxDeg > 0:
        amb = [a for a in amb if a.deg-2 <= MaxDeg]

    return amb
############################################################################
def generateOverlaps(a,b):
    v = a[0]
    w = b[0]
    return [Overlap(v + w[k:],v[:-k+1],w[k-1:],a[1],b[1]) for k in range(2,min(len(v),len(w))-1) if v.endswith(w[:k])]
############################################################################
def generateInclusions(a,b):
    v = a[0]
    w = b[0]
    k = 0
    if len(w) <= len(v) and a[1] != b[1]:
        while True:
            k = v.find(w, k)+1
            if k > 0:
                yield Inclusion(v,v[:k],v[k+len(w)-2:],a[1],b[1])
            else:
                break
############################################################################
def DeleteRedundant(amb,lt,Info=True):
    start = time()
    overlaps = [a for a in amb if type(a) == Overlap]
    incls = [a for a in amb if type(a) == Inclusion]
    
    for V,k in lt:
        overlaps = [a for a in overlaps if not (k < a.min and V in a.ABC)]
        incls = [a for a in incls if not (k < a.j and _incl_test_(a,V))]
    
    result = overlaps + incls
    if Info:
        print("Removed " + str(len(amb) - len(result)) + " ambiguities in %.5f" % (time()-start))
    return result

def _incl_test_(a,V):
    """
    V|A or V|B or V|C or B|V|ABC but V != ABC
    """
    if V in a.A or V in a.C:
        return True
    B = a.ABC[len(a.A)-1 : -len(a.C)+1]
    if V in B or (B in V and V in a.ABC[1:-1]):
        return True
    else:
        return False
