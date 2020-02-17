# coding: utf-8
"""
F4
================

Module that implements the F4 algorithm in the free algebra that
allows to trace cofactors.

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
from collections import defaultdict, deque

from sage.all import *
from sage.matrix.matrix_rational_sparse import Matrix_rational_sparse

from .Ambiguities import Overlap,Inclusion,GenerateAmbiguities,DeleteRedundant,generateInclusions
from .NCMonomial import NCMonomial
from .NCPoly import NCPoly
from .NormalForms import __interreduce__

############################################################################
# S-polynomials
############################################################################
class CritPair:
    def __init__(self,amb,fi,fj):
        self.deg = amb.deg
        one = NCMonomial.one()
        mon_a = NCMonomial(1,amb.A,intern=True)
        mon_c = NCMonomial(1,amb.C,intern=True)
        
        if isinstance(amb,Inclusion):
            f2 = fj.lrmul(amb.A,amb.C)
            if fi == f2:
                self.deg = -1
            else:
                self.f = fi
                self.g = f2
                self.triple_f = (one,amb.i,one)
                self.triple_g = (mon_a,amb.j,mon_c)
        else:
            f1 = fi.rmul(amb.C)
            f2 = fj.lmul(amb.A)
            if f1 == f2:
                self.deg = -1
            else:
                self.f = f1
                self.g = f2
                self.triple_f = (one,amb.i,mon_c)
                self.triple_g = (mon_a,amb.j,one)
############################################################################
def CheckResolvability(G,oldlength=0,Criterion=True,MaxDeg=-1,Info=True):
    
    words_new = [(g.lt.mon,i+oldlength) for (i,g) in enumerate(G[oldlength:])]
    amb = []
    
    for ((m,i),(j,g)) in itertools.product(words_new,enumerate(G[:oldlength])):
        if not g.redundant and m in g.lt.mon:
            g.redundant = True
            incl = next(generateInclusions((g.lt.mon,j),(m,i)))
            amb.append(incl)
    
    words_old = [(g.lt.mon,i) for (i,g) in enumerate(G[:oldlength]) if not g.redundant]
    start = time()
    amb += GenerateAmbiguities(words_old,words_new,MaxDeg=MaxDeg)
    if Info:
        print(str(len(amb)) + " ambiguities in total (computation took %.5f)" % (time()-start))
    
    if Criterion:
        amb = DeleteRedundant(amb,words_old + words_new,Info=Info)

    critPairs = [CritPair(a,G[a.i],G[a.j]) for a in amb]
    critPairs = [pair for pair in critPairs if pair.deg != -1]
    critPairs.sort(key=lambda p: p.deg)

    if Info:
        print(str(len(critPairs)) + " critical pairs were generated.")

    return critPairs
############################################################################
# Groebner basis computation
############################################################################
def Groebner(cofactors,F,MaxIter=10,Ignore=0,MaxDeg=-1,Criterion=True,Info=False,IterCount=0):
    intern = isinstance(F[0],NCPoly)
    if intern:
        G = F
    else:
        G = [NCPoly(f) for f in F]
        for i,g in enumerate(G):
            g.makeMonic(i)

    oldlength = len(G)
    count = 0
    if Info:
        print("G has " + str(len(G)) + " elements in the beginning.")

    critPairs = CheckResolvability(G,Ignore,Criterion=Criterion,MaxDeg=MaxDeg,Info=Info)

    while count < MaxIter:
        start = time()
        while len(critPairs) > 0:
            critPairs.sort(key = lambda p : p.deg)
            d = critPairs[0].deg
            P = [pair for pair in critPairs if pair.deg == d]
            critPairs = critPairs[len(P):]
            PPlus = Reduction(P,G)
            G += PPlus
            print(str(len(critPairs)).ljust(10))
            sys.stdout.write("\033[F")
   
        end = time()
        if Info:
            print("Reduction took: %.5f" % (end-start))
        if len(G) == oldlength:
            if Info:
                print("All critical pairs could be reduced to 0.\n")
            break

        if Info:
            print("Iteration " + str(count+IterCount+1) + " finished. G has now " + str(len(G)) + " elements.\n")
        
        count += 1
        if count < MaxIter:
            P = G[oldlength:]
            __prepForInterreduction__(P,oldlength)
            G = G[:oldlength] + __interreduce__(P)
            critPairs = CheckResolvability(G,oldlength,Criterion=Criterion,MaxDeg=MaxDeg,Info=Info)
            oldlength = len(G)

    if intern:
        return G
    else:
        P = F[0].parent()
        cofactors[:] = [c for i,c in enumerate(__rewriteCofactors__(G,F,Info=Info)) if not G[i].redundant]
        return [g.toNormal(P) for g in G if not g.redundant]
############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
def __prepForInterreduction__(P,oldlength):
    for p in P:
        bad_triples = ((a,i,b) for (a,i,b) in p.cofactors if i > oldlength)
        cofactors = [(a,i,b) for (a,i,b) in p.cofactors if i <= oldlength]
        for (a,i,b) in bad_triples:
            cofactors += [(a*l,j,r*b) for (l,j,r) in P[i-oldlength].cofactors]
        p.cofactors = cofactors
############################################################################
def __rewriteCofactors__(G,F,Info=False):
    if Info:
        print("Rewriting the cofactors has started.")
    start = time()
    P = F[0].parent()
    N = len(F)

    cofactors = [[(a.toNormal(P),j,b.toNormal(P)) for (a,j,b) in g.cofactors] for g in G]
    SimplifyCofactors(cofactors)
    
    for i,c in enumerate(cofactors[:N]):
        cofactors[i] = [(a,F[j],b) for (a,j,b) in c]
    
    for idx,c in enumerate(cofactors[N:]):
        j = idx + N
        cofactors_temp = []
        for (a,i,b) in c:
            cofactors_temp += [(a*l,f,r*b) for (l,f,r) in cofactors[i]]
        cofactors[j] = cofactors_temp
    end = time()
    if Info:
        print("Rewriting the cofactors took in total %.5f\n" %(end-start))
    return cofactors
############################################################################
def Reduction(P,G):
    F = [f for pair in P for f in (pair.f,pair.g)]
    pivot_rows,pivot_cols,columns,cofactors_G = SymbolicPreprocessing(F,G)

    #split ciritcal polynomials in pivot and non-pivot rows
    rest_rows = [pair.g for pair in P]
    for pair in P:
        f = pair.f
        if f.lt.mon in pivot_cols:
            rest_rows.append(f)
        else:
            pivot_cols.add(f.lt.mon)
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
    pivot_rows.sort(key=lambda f: f.lt,reverse=True)
    rows = pivot_rows + rest_rows

    (A,B,C,D) = getMatrices(rows,columns,n,m)

    A_inv = A.rref()[:,n:]
    CA_inv = C * A_inv
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
    
    FPlus = [NCPoly(f[0],f[1:],[]) for f in FPlus]
    
    # take care of cofactors
    for pair in P:
        cofactors_G[pair.f] = [pair.triple_f]
        cofactors_G[pair.g] = [pair.triple_g]

    cofactorsF = [cofactors_G[f] for f in rows]
    T1 = -T * CA_inv
    for (i,j) in T1[:rank,:].nonzero_positions():
        coeff = T1[i,j]
        FPlus[i].cofactors += [multiplyTriple(t,coeff) for t in cofactorsF[j]]

    for (i,j) in T[:rank,:].nonzero_positions():
        coeff = T[i,j]
        FPlus[i].cofactors += [multiplyTriple(t,coeff) for t in cofactorsF[j+n]]

    FPlus.sort(key=lambda f : f.lt)

    return FPlus
############################################################################
def SymbolicPreprocessing(F,G):
    columns = {f.lt.mon for f in F}
    T = {m.mon for f in F for m in f.tail} - columns
    lt = [(i,g.lt.mon) for i,g in enumerate(G) if not g.redundant]
    G_prime = []
    cofactors_G = dict()
    G_mons = []
    while len(T) > 0:
        t = T.pop()
        columns.add(t)
        (i,m) = next(((i,m) for i,m in lt if m in t), (None,None))
        if m:
            g = G[i]
            factors = t.split(m,1)
            a = factors[0] + '*'
            b = '*' + factors[1]
            agb = g.lrmul(a,b)
            G_prime.append(agb)
            G_mons.append(agb.lt.mon)
            cofactors_G[agb] = [(NCMonomial(1,a,intern=True),i,NCMonomial(1,b,intern=True))]
            T.update({m.mon for m in agb.tail if m.mon not in columns})
    return G_prime,set(G_mons),columns,cofactors_G
############################################################################
def getMatrices(rows,columns,n,m):
    entries = {(i,i+n):1 for i in range(n)}
    for i,f in enumerate(rows):
        entries[(i,columns[f.lt.mon])] = f.lt.coeff
        entries.update({(i,columns[m.mon]):m.coeff for m in f.tail})
    l = len(columns) + n
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
def SimplifyCofactors(cofactors):
    for idx,c in enumerate(cofactors):
        already_checked = defaultdict(deque)
        to_remove = []
        for i,(a,f,b) in enumerate(c):
            is_in = already_checked[(-a,f,b)]
            if is_in:
                j = is_in.pop()
                to_remove += [i,j]
            else:
                already_checked[(a,f,b)].append(i)
        to_remove = set(to_remove)
        cofactors[idx] = [t for j,t in enumerate(c) if j not in to_remove]
