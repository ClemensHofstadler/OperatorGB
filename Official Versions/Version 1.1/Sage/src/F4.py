# coding: utf-8
"""
F4
================

Module that implements the F4 algorithm in the free algebra that
allows to trace cofactors.

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
from time import time

from sage.all import *
from sage.matrix.matrix_rational_sparse import Matrix_rational_sparse

from Ambiguities import Overlap,Inclusion,GenerateAmbiguities,DeleteRedundant
from NCMonomial import NCMonomial
from NCPoly import NCPoly

############################################################################
# S-polynomials
############################################################################
class CritPair:
    def __init__(self,amb,fi,fj):
        one = NCMonomial(1,'*',intern=True)
        mon_a = NCMonomial(1,amb.A,intern=True)
        mon_c = NCMonomial(1,amb.C,intern=True)
        self.deg = amb.deg
        if isinstance(amb,Inclusion):
            f = fj.lrmul(amb.A,amb.C)
            if fi == f:
                self.deg = -1
            else:
                self.f = fi
                self.g = f
                self.triple_f = (one,fi,one)
                self.triple_g = (mon_a,fj,mon_c)
        else:
            f1 = fi.rmul(amb.C)
            f2 = fj.lmul(amb.A)
            if f1 == f2:
                self.deg = -1
            else:
                self.f = f1
                self.g = f2
                self.triple_f = (one,fi,mon_c)
                self.triple_g = (mon_a,fj,one)
############################################################################
def CheckResolvability(G,oldlength=0,Criterion=True,MaxDeg=-1,Info=True):
    words = [(f.lt.mon,i) for (i,f) in enumerate(G)]
    start = time()
    amb = GenerateAmbiguities(words[:oldlength],words[oldlength:],MaxDeg=MaxDeg)
    if Info:
        print str(len(amb)) + " ambiguities in total (computation took %.5f)" % (time()-start)
    
    if Criterion:
        amb = DeleteRedundant(amb,[f.lt.mon for f in G],Info=Info)

    critPairs = [CritPair(a,G[a.i],G[a.j]) for a in amb]
    critPairs = [pair for pair in critPairs if pair.deg != -1]
    critPairs.sort(key=lambda p: p.deg)

    if Info:
        print str(len(critPairs)) + " critical pairs were generated."

    return critPairs
############################################################################
# Groebner basis computation
############################################################################
def Groebner(cofactors,ideal,MaxIter=10,N=-1,Ignore=0,MaxDeg=-1,Criterion=True,Info=True,OutputProd=False,IterCount=0):
    if not isinstance(ideal[0],NCPoly):
        G = map(NCPoly,ideal)
    else:
        G = ideal

    for f in G:
        f.makeMonic()
    
    oldlength = len(G)
    count = 0
    if Info:
        print "G has " + str(len(G)) + " elements in the beginning."

    critPairs = CheckResolvability(G,Ignore,Criterion=Criterion,MaxDeg=MaxDeg,Info=Info)

    while count < MaxIter:
        start = time()
        while len(critPairs) > 0:
            d = critPairs[0].deg
            F = [pair for pair in critPairs if pair.deg == d]
            critPairs = critPairs[len(F):]
            FPlus = Reduction(F,G)
            F = [f for (f,c) in FPlus]
            G += F
            cofactors += FPlus
            print str(len(critPairs)).ljust(10)
            sys.stdout.write("\033[F")
   
        end = time()
        if Info:
            print "Reduction took: %.5f" % (end-start)
        if len(G) == oldlength:
            print "All critical pairs could be reduced to 0.\n"
            break
        
        count += 1
        if Info:
            print "Iteration " + str(count+IterCount) + " finished. G has now " + str(len(G)) + " elements.\n"
        if count < MaxIter:
            critPairs = CheckResolvability(G,oldlength,Criterion=Criterion,MaxDeg=MaxDeg,Info=Info)
            oldlength = len(G)
    
    if OutputProd:
        return G
    else:
        start = time()
        if Info:
            print "Rewriting the cofactors has started."
        __rewriteCofactors__(cofactors,ideal)
        end = time()
        if Info:
            print "Rewriting the cofactors took in total %.5f\n" %(end-start)
        return [f.toNormal() for f in G]
############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
def Reduction(L,G):

    F = [f for pair in L for f in (pair.f,pair.g)]
    (pivot_rows,pivot_cols,columns,cofactors_G) = SymbolicPreprocessing(F,G)
    
    #split ciritcal polynomials in pivot and non-pivot rows
    rest_rows = [pair.g for pair in L]
    for pair in L:
        f = pair.f
        if f.lt.mon in pivot_cols:
            rest_rows.append(f)
        else:
            pivot_cols.append(f.lt.mon)
            pivot_rows.append(f)
    #seperate pivot and non-pivot columns
    rest_cols = [NCMonomial(1,m,intern=True) for m in columns if m not in pivot_cols]
    pivot_cols = [NCMonomial(1,m,intern=True) for m in pivot_cols]
    pivot_cols.sort(reverse=True)
    rest_cols.sort(reverse=True)
    n = len(pivot_rows)
    m = len(rest_rows)

    columns = {m.mon:i for (i,m) in enumerate(pivot_cols)}
    columns.update({m.mon:i+2*n for (i,m) in enumerate(rest_cols)})
    rows = pivot_rows + rest_rows

    (A,B,C,D) = getMatrices(rows,columns,n,m)
    
    A_inv = A.rref()[:,n:]
    CA_inv = C*A_inv
    
    #bring D - CA^{-1}B in RRef
    diff(D,CA_inv*B)
    M = D.rref()
    #get transformation matrix
    T = M[:,-m:]
    M = M[:,:-m]

    #compute polynomials
    pos = M.nonzero_positions()
    if len(pos) == 0:
        return []
    rank = pos[-1][0]+1
    FPlus = [[] for i in range(rank)]
    for (i,j) in pos:
        FPlus[i].append(NCMonomial(M[i,j],rest_cols[j].mon,intern=True))

    #compute cofactors
    cofactors = [[] for i in range(rank)]
    for pair in L:
        cofactors_G[pair.f] = pair.triple_f
        cofactors_G[pair.g] = pair.triple_g
    cofactors_F = [cofactors_G[f] for f in rows]
    T1 = -T*CA_inv
    for (i,j) in T1[:rank,:].nonzero_positions():
        triple = copyTriple(cofactors_F[j])
        triple[0].coeff *= T1[i,j]
        cofactors[i].append(triple)
    for (i,j) in T[:rank,:].nonzero_positions():
        triple = copyTriple(cofactors_F[j+n])
        triple[0].coeff *= T[i,j]
        cofactors[i].append(triple)

    FPlus = [(NCPoly(f[0],f[1:]),c) for (f,c) in zip(FPlus,cofactors)]

    return FPlus
############################################################################
def SymbolicPreprocessing(F,G):
    columns = {f.lt.mon for f in F}
    T = {m.mon for f in F for m in f.tail}
    lt = [(g.lt.mon,g) for g in G]
    G_prime = []
    cofactors_G = dict()
    G_mons = []
    while len(T) > 0:
        t = T.pop()
        columns.add(t)
        for (m,g) in lt:
            if m in t:
                factors = t.split(m,1)
                l = factors[0] + '*'
                r = '*' + factors[1]
                f = g.lmul(l).rmul(r)
                G_prime.append(f)
                G_mons.append(f.lt.mon)
                cofactors_G[f] = (NCMonomial(1,l,intern=True),g,NCMonomial(1,r,intern=True))
                T.update({m.mon for m in f.tail if m.mon not in columns})
                break
    return (G_prime,G_mons,columns,cofactors_G)
############################################################################
def getMatrices(rows,columns,n,m):
    entries = {(i,i+n):1 for i in range(n)}
    for i,f in enumerate(rows):
        entries[(i,columns[f.lt.mon])] = f.lt.coeff
        entries.update({(i,columns[m.mon]):m.coeff for m in f.tail})
    l = len(columns)+n
    entries.update({(i+n,i+l):1 for i in range(m)})
    
    M = matrix(QQ,entries,sparse=True)

    A = M[:n,:2*n]
    B = M[:n,2*n:len(columns)+n]
    C = M[n:,:n]
    D = M[n:,2*n:]

    return (A,B,C,D)
############################################################################
def diff(D,M):
    for (i,j) in M.nonzero_positions():
        D.add_to_entry(i,j,-M[i,j])
############################################################################
# Rewriting cofactors
############################################################################
def __rewriteCofactors__(cofactors,I):
    start = time()
    cofactors[:] = [(f.toNormal(),[tripleToNormal(t) for t in c]) for (f,c) in cofactors]
    end = time()
    cofactors_done = {cofactors[0][0]:cofactors[0][1]}
    print "convert data structure %.5f" %(end-start)
    
    for (i,(f,c)) in enumerate(cofactors):
        if i == 0:
            continue
        cofactors_temp = []
        for (a,g,b) in c:
            if g in cofactors_done:
                cofactors_g = cofactors_done[g]
                cofactors_temp += [(a*l,h,r*b) for (l,h,r) in cofactors_g]
            else:
                cofactors_temp.append((a,g,b))
        
        cofactors[i] = (f,cofactors_temp)
        cofactors_done[f] = cofactors_temp
############################################################################
# additional stuff
############################################################################
def tripleToNormal(t):
    return (t[0].toNormal(),t[1].toNormal(),t[2].toNormal())
############################################################################
def copyTriple(t):
    return (t[0].copy(),t[1].copy(),t[2].copy())
############################################################################
def multiplyTriple(triple,c):
    f = copyTriple(triple)
    f[0].coeff *= c
    return f
