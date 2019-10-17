# coding: utf-8
"""
Ambiguites
================

Module that provides functionality to compute ambiguities of words

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

import itertools
from time import time

from MonomialOrder import SortedQ

############################################################################
# ambiguities
############################################################################
class Overlap:
    """
    fiC - Afj
    """
    def __init__(self,ABC,C,A,i,j):
        self.ABC = ABC
        self.C = C
        self.A = A
        self.i = i
        self.j = j
        self.deg = len(ABC.split('*'))
        self.min = min(i,j)
        self.max = max(i,j)

    def __repr__(self):
        return "Overlap(" + str(self.ABC) + ", " + str(self.C) + ", " + str(self.A) + ", (" + str(self.i) + ", " + str(self.j) + "))"
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
        self.max = max(i,j)

    def __repr__(self):
        return "Inclusion(" + str(self.ABC) + ", " + str(self.A) + ", " + str(self.C) + ", (" + str(self.i) + ", " + str(self.j) + "))"
############################################################################
def GenerateAmbiguities(words,newPart = None,MaxDeg = -1):
    if newPart != None:
        return __generateAmbiguities__(words,newPart,MaxDeg)
    
    amb = [generateOverlaps(v,w) for (v,w) in itertools.product(words,repeat=2)]
    amb += [generateInclusions(v,w) for (v,w) in itertools.product(words,repeat=2)]

    amb = flatten(amb)
   
    if MaxDeg > 0:
        amb = [a for a in amb if a.deg-2 <= MaxDeg]

    return amb
############################################################################
def __generateAmbiguities__(oldPart,newPart,MaxDeg):
    
    amb = [generateOverlaps(v,w) + generateOverlaps(w,v) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [generateInclusions(v,w) + generateInclusions(w,v) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [generateOverlaps(v,w) for (v,w) in itertools.product(newPart,repeat=2)]
    amb += [generateInclusions(v,w) for (v,w) in itertools.product(newPart,repeat=2)]
    
    amb = flatten(amb)
    
    if MaxDeg > 0:
        amb = [a for a in amb if a.deg-2 <= MaxDeg]

    return amb
############################################################################
def generateOverlaps(a,b):
    v = a[0]
    w = b[0]
    overlaps = []
    for k in range(2,min(len(v),len(w))-1):
        if v.endswith(w[:k]):
            overlaps.append(Overlap(v + w[k:],w[k-1:],v[:-k+1],a[1],b[1]))
    return overlaps
############################################################################
def generateInclusions(a,b):
    v = a[0]
    w = b[0]
    k = 0
    inclusions = []
    if len(w) < len(v):
        while True:
            k = v.find(w, k)+1
            if k > 0:
                inclusions.append(Inclusion(v,v[:k],v[k+len(w)-2:],a[1],b[1]))
            else:
                break
    return inclusions
############################################################################
def DeleteRedundant(amb,lt,Info=True):
    start = time()
    amb.sort(key=lambda a:a.deg,reverse=True)
    indices = [k for a in amb for k in (a.i,a.j)]
    i = 0
    for k in iter(set(indices)):
        selected = [a for a in amb if a.max == k]
        while len(selected) > 0:
            f = selected.pop()
            j = f.min
            ABC = f.ABC
            #Chain criterion: we can remove f if there is V in lt such that V|ABC and index(V) < f.min
            #So, we have to keep f if V not in ABC for all V in lt[:f.min]
            if all(V not in ABC for V in lt[:j]):
                yield f
                i += 1
            selected = [a for a in selected if ABC not in a.ABC or (a.min <= j and ABC == a.ABC and (a.min != j or SortedQ(a.A,f.A)))]
    if Info:
        print "Removed " + str(len(amb) - i) + " ambiguities in %.5f" % (time()-start)
############################################################################
# additional stuff
############################################################################
def flatten(l):
    return [item for sublist in l for item in sublist]
