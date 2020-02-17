# coding: utf-8
"""
OperatorGB
================

Module to compute normal forms of noncommutative polynomials

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
from random import choice

from sage.all import *
from sage.matrix.matrix_rational_sparse import Matrix_rational_sparse

from .Ambiguities import flatten
from .NCMonomial import NCMonomial
from .NCPoly import NCPoly

############################################################################
# Normal form computation
############################################################################
def ReducedForm(cofactorsF,G_input,f_input):
    G = [NCPoly(g) for g in G_input]
    for i,g in enumerate(G):
        g.makeMonic(i)
    f = NCPoly(f_input)
    f.cofactors = [(NCMonomial.one(),len(G),NCMonomial.one())]

    normal_form = __reducedForm__(G+[f],len(G),intern=False)

    cofactorsF[:] = [tripleToNormal((a.__mul__(-1),G[i],b)) for (a,i,b) in normal_form.cofactors[1:]]
    return normal_form.toNormal()
############################################################################
def Rewrite(cofactorsF, cofactorsG):
    G = list(map(MultiplyOut,cofactorsG))
    result = []
    for (a,g,b) in cofactorsF:
        i = G.index(g)
        result += [(a*l,f,r*b) for (l,f,r) in cofactorsG[i]]
    return result

############################################################################
def __reducedForm__(G,i,intern=True):
    f = G[i]
    F = [f]
    columns,cofactors = __symbolicPreprocessingRed__(F,G,i)
    cofactors = [(NCMonomial.one(),i,NCMonomial.one())] + cofactors
    columns = [NCMonomial(1,c,intern=True) for c in columns]
    columns.sort(reverse=True)

    M = SetUpMatrix(F,columns).rref()
    T = M[:,len(columns):] #transformation matrix
    RREF = M[:,:len(columns)]

    row_idx = T.nonzero_positions_in_column(0)[-1]

    # reduction to zero
    if len(RREF.nonzero_positions_in_row(row_idx)) == 0:
        if intern:
            return 0
        else:
            normal_form = NCPoly.zero()
    else:
        monomials = []
        for j in RREF.nonzero_positions_in_row(row_idx):
            monomials.append(NCMonomial(RREF[row_idx,j],columns[j].mon,intern=True))
        normal_form = NCPoly(monomials[0],monomials[1:],[])

    # no reduction
    if intern and normal_form == f:
        return None
    
    # some reduction - take care of cofactors
    coeff = T[row_idx,0]
    for j in T.nonzero_positions_in_row(row_idx):
        (a,idx,b) = cofactors[j]
        coeff = T[row_idx,j]
        normal_form.cofactors += [(a*l*coeff,k,r*b) for (l,k,r) in G[idx].cofactors]
    return normal_form

############################################################################
def __symbolicPreprocessingRed__(F,G,idx,shuffle=False):
    T = {m.mon for m in F[0].tail}
    T.add(F[0].lt.mon)
    lt = [(i,g.lt.mon) for i,g in enumerate(G) if i != idx and g is not None]
    done = set()
    cofactors = []
    while len(T) > 0:
        t = T.pop()
        done.add(t)
        if shuffle:
            choices = [(i,m) for i,m in lt if m in t]
            i,m = choice(choices) if choices else (None,None)
        else:
            i,m = next(((i,m) for i,m in lt if m in t), (None,None))
        if m:
            g = G[i]
            factors = t.split(m,1)
            a = factors[0] + '*'
            b = '*' + factors[1]
            agb = g.lrmul(a,b)
            cofactors.append([NCMonomial(1,a,intern=True),i,NCMonomial(1,b,intern=True)])
            F.append(agb)
            T.update({m.mon for m in agb.tail if m.mon not in done})
    return done, cofactors
############################################################################
def __idealMembership__(G,f,shuffle=False):
    F = [f]
    columns,cofactors = __symbolicPreprocessingRed__(F,G,-1,shuffle=shuffle)
    columns = [NCMonomial(1,c,intern=True) for c in columns]
    columns.sort(reverse=True)
    M = SetUpMatrix(F[1:],columns)[:,:len(columns)]
    b = vector(SetUpMatrix(F[:1],columns).list()[:-1])
    try:
        coeffs = M.solve_left(b)
    except:
        return False

    cofactors = [multiplyTriple(t,c) for c,t in zip(coeffs,cofactors) if c != 0]

    return cofactors
############################################################################
def SetUpMatrix(F,columns):
    entries = {}
    l = len(columns)
    cols = {c.mon:i for i,c in enumerate(columns)}
    for i,f in enumerate(F):
        entries[(i,cols[f.lt.mon])] = f.lt.coeff
        entries[(i,i+l)] = 1
        entries.update({(i,cols[m.mon]):m.coeff for m in f.tail})

    return matrix(QQ,entries,sparse=True)
############################################################################
# additional stuff
############################################################################
def tripleToNormal(t):
    return (t[0].toNormal(),t[1].toNormal(),t[2].toNormal())
############################################################################
def copyTriple(t):
    return (t[0].copy(),t[1],t[2].copy())
############################################################################
def multiplyTriple(triple,c):
    t = copyTriple(triple)
    t[0].coeff *= c
    return t
############################################################################
# Interreduction
############################################################################
def __interreduce__(G):
    i = 0
    s = len(G)
    while i < s:
        if G[i] is None:
            i += 1
            continue
    
        normal_form = __reducedForm__(G,i)
        # normal form = 0
        if type(normal_form) is int:
            G[i] = None
        # normal form != 0
        elif normal_form:
            G[i] = normal_form
            i = 0
        # normal form = g_i
        else:
            i += 1

    return [g for g in G if g]
############################################################################
# Checking correctness
############################################################################
def MultiplyOut(cofactors):
    return sum(map(prod,cofactors))
############################################################################
def AllInI(cofactors,I):
    return all([g in I for (l,g,r) in cofactors])
############################################################################
def CheckCofactors(cofactors,G,I):
    print("cofactors = G:")
    print(all(MultiplyOut(c) == G[i] for i,c in enumerate(cofactors)))

    print("linear combination only of elements in I:")
    print(all(AllInI(c,I) for c in cofactors))
############################################################################
def CheckCertificate(certificate,claim,I):
    print("certificate = claim:")
    print(MultiplyOut(certificate) == claim)

    print("certificate consists only of elements in I:")
    print(all(g in I for (l,g,r) in certificate))
############################################################################
def _checkInterreduction_(I,G):
    P = G[0].parent()
    cofactors = []
    for f in I:
        cofactors.append([(a.toNormal(P),G[j],b.toNormal(P)) for (a,j,b) in f.cofactors])
    print(all(MultiplyOut(c) == f.toNormal(P) for c,f in zip(cofactors,I)))
