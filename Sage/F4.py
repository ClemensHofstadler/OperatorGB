# coding: utf-8
"""
F4
================

Package to compute noncommutative Gröbner bases in the free Algebra over QQ
together with tracing of cofactors.

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

from sage.matrix.matrix_rational_sparse import Matrix_rational_sparse

import itertools
from time import time
from random import sample

from sage.all import *

############################################################################
# global variables
############################################################################
WordOrder = []
Knowns = []
Unknowns = []
SortedQ = None
t1 = t2 = t3 = t4 = t5 = t6 = 0
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
        F = FreeAlgebra(QQ,len(WordOrder),WordOrder)
        return F(s)
############################################################################
# NCPoly
############################################################################
class NCPoly:
    def __init__(self,f,*args):
        if args:
            self.lt = f
            self.tail = args[0]
        else:
            d = f.monomial_coefficients()
            if d:
                monomials = [NCMonomial(d[key],key) for key in d.keys()]
                monomials = sorted(monomials)
            else:
                monomials = [NCMonomial(0,'*',intern=True)]
            self.lt = monomials[-1]
            self.tail = monomials[:-1]
############################################################################
    def copy(self):
        return NCPoly(self.lt.copy(),[mon.copy() for mon in self.tail])
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
        if len(self.tail) != len(other.tail):
            return False
        return self.lt == other.lt and all([m1 == m2 for (m1,m2) in zip(self.tail,other.tail)])
############################################################################
    def __hash__(self):
        return hash((self.lt,) + tuple(self.tail))
############################################################################
    def makeMonic(self):
        lc = self.lt.coeff
        if lc != 1:
            self.lt.coeff = 1
            for m in self.tail:
                m.coeff /= lc
############################################################################
    def lmul(self,m):
        f = self.copy()
        f.lt.mon = m[:-1] + f.lt.mon
        for mon in f.tail:
            mon.mon = m[:-1] + mon.mon
        return f
############################################################################
    def rmul(self,m):
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
    def toNormal(self):
        return self.lt.toNormal() + sum([_.toNormal() for _ in self.tail])
############################################################################
# set up the ring
############################################################################
def SetUpRing(vars,unknowns = None,Info=True):
    if unknowns != None:
        return __SetUpRing__(vars,unknowns,Info=Info)
    
    global WordOrder, SortedQ
    WordOrder = [str(v) for v in vars]
    SortedQ = DegLex
    if Info:
        s = ""
        for v in WordOrder:
            s += v + " < "
        print s[:-3]
############################################################################
def __SetUpRing__(knowns,unknowns,Info=True):
    global WordOrder, Knowns, Unknowns, SortedQ
    Knowns = [str(v) for v in knowns]
    Unknowns = [str(v) for v in unknowns]
    WordOrder = Knowns + Unknowns
    SortedQ = MultiLex
    if Info:
        s = ""
        for v in Knowns:
            s += v + " < "
        s = s[:-1]
        s += "< "
        for v in Unknowns:
            s += v + " < "
        print s[:-3]
############################################################################
# monomial order
############################################################################
def DegLex(a,b):
    a_vars = a.split('*')
    b_vars = b.split('*')
    la = len(a_vars)
    lb = len(b_vars)
    
    if la == lb:
        i = 0
        while i < la and a_vars[i] == b_vars[i]:
            i+=1
        if i == la:
            return True
        else:
            return WordOrder.index(a_vars[i]) < WordOrder.index(b_vars[i])

    else:
        return la < lb
############################################################################
def MultiLex(a,b):
    a_vars = a.split('*')
    b_vars = b.split('*')
    V1a = len([v for v in a_vars if v in Unknowns])
    V2a = len([v for v in a_vars if v in Knowns])
    V1b = len([v for v in b_vars if v in Unknowns])
    V2b = len([v for v in b_vars if v in Knowns])

    if V1a < V1b or (V1a == V1b and V2a < V2b):
        return True
    if V1a > V1b or (V1a == V1b and V2a > V2b):
        return False

    return DegLex(a,b)
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

    def __repr__(self):
        return "Inclusion(" + str(self.ABC) + ", " + str(self.A) + ", " + str(self.C) + ", (" + str(self.i) + ", " + str(self.j) + "))"
############################################################################
def GenerateAmbiguities(words,newPart = None,MaxDeg = -1):
    if newPart != None:
        return __GenerateAmbiguities__(words,newPart,MaxDeg)
    
    amb = [generateOverlaps(v,w) for (v,w) in itertools.product(words,repeat=2)]
    amb += [generateInclusions(v,w) for (v,w) in itertools.product(words,repeat=2)]

    amb = flatten(amb)
   
    if MaxDeg > 0:
        amb = [a for a in amb if a.deg-2 <= MaxDeg]

    return amb
############################################################################
def __GenerateAmbiguities__(oldPart,newPart,MaxDeg):
    
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
def DeleteRedundant(amb,Info=True):
    start = time()
    amb.sort(key=lambda a:a.deg,reverse=True)
    indices = [k for a in amb for k in (a.i,a.j)]
    i=0
    for k in iter(set(indices)):
        selected = [a for a in amb if max(a.i,a.j)==k]
        while len(selected) > 0:
            f = selected.pop()
            yield f
            i+=1
            selected = [a for a in selected if not f.ABC in a.ABC]
    end = time()
    if Info:
        print("Removed " + str(len(amb)-i) + " ambiguities in %.5f" % (end-start))
############################################################################
# S-polies
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
    end = time()
    if Info:
        print(str(len(amb)) + " ambiguities in total (computation took %.5f)" % (end-start) )
    
    if Criterion:
        amb = DeleteRedundant(amb,Info=Info)

    critPairs = [CritPair(a,G[a.i],G[a.j]) for a in amb]
    critPairs = [pair for pair in critPairs if pair.deg != -1]
    critPairs.sort(key=lambda p: p.deg)

    if Info:
        print(str(len(critPairs)) + " critical pairs were generated.")

    return critPairs
############################################################################
# F4
############################################################################
def F4(cofactors,ideal,MaxIter=10,N=-1,Ignore=0,MaxDeg=-1,Criterion=True,Info=True,OutputProd=False,IterCount=0):
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

    global t1,t2,t3,t4,t5,t6
    
    start = time()
    F = [f for pair in L for f in (pair.f,pair.g)]
    (pivot_rows,pivot_cols,columns,cofactors_G) = SymbolicPreprocessing(F,G)
    t1 += time() -start
    
    #split ciritcal polynomials in pivot and non-pivot rows
    start = time()
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
    t2+= time() - start

    start = time()
    (A,B,C,D) = getMatrices(rows,columns,n,m)
    t3 += time()-start
    
    start = time()
    A_inv = A.rref()[:,n:]
    CA_inv = C*A_inv
    t4 += time() - start
    
    start = time()
    #bring D - CA^{-1}B in RRef
    diff(D,CA_inv*B)
    M = D.rref()
   #get transformation matrix
    T = M[:,-m:]
    M = M[:,:-m]
    t5 += time()-start
    
    
    start = time()
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
    t6 += time() -start

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
def multiplyTriple(triple,c):
    f = copyTriple(triple)
    f[0].coeff *= c
    return f
############################################################################
# Normal form computation
############################################################################
def ReducedForm(G_Input,f,InputProd=False):
    if not InputProd:
        G = map(NCPoly,G_Input)
        p = NCPoly(f)
    else:
        G = G_Input
        p = f

    cofactors = []
    (F,columns,cofactorsF) = __symbolicPreprocessingRed__(p,G)
    columns = sorted(columns,reverse=True)
    cofactorsF = [(NCMonomial(1,'*',intern=True),p,NCMonomial(1,'*'))] + cofactorsF
    F = [p] + F
    M = SetUpMatrix(F,columns)
    M = M.rref()
    T = M[:,len(columns):] #transformation matrix
    R = M[:,:len(columns)]
    
    pos = T.nonzero_positions_in_column(0)[-1]
    coeffs = T[pos,:]
    coeff_f = coeffs[0,0]
    coeffs = [c/coeff_f for c in coeffs.list()]

    #reduction to 0
    if len(R.nonzero_positions_in_row(pos)) == 0:
        normal_form = 0
        if InputProd:
            cofactors = [(NCMonomial(-c*l.coeff,l.mon,intern=True),g,r) for (c,(l,g,r)) in zip(coeffs[1:],cofactorsF[1:]) if c != 0]
        else:
            cofactors = [(-c*l.toNormal(),g.toNormal(),r.toNormal()) for (c,(l,g,r)) in zip(coeffs[1:],cofactorsF[1:]) if c != 0]

    #reduction to normal form != 0
    else:
        (normal_form,cofactors) = Matrix2Poly(M[pos,:],columns,[],cofactorsF)[0]
        normal_form = normal_form.toNormal()/coeff_f
        if InputProd:
            normal_form = NCPoly(normal_form)
            cofactors = [(NCMonomial(-l.coeff/coeff_f,l.mon,intern=True),g,r) for (l,g,r) in cofactors[1:]]
        else:
            cofactors = [(-l.toNormal()/coeff_f,g.toNormal(),r.toNormal()) for (l,g,r) in cofactors[1:]]

    return (normal_form,cofactors)
############################################################################
def IdealMembership(G_Input,f,InputProd=False):
    if not InputProd:
        G = map(NCPoly,G_Input)
    else:
        G = G_Input

    p = NCPoly(f)
    cofactors = []
    (F,columns,cofactors) = __symbolicPreprocessingRed__(p,G)
    columns = sorted(columns,reverse=True)
    M = SetUpMatrix(F,columns)[:,:len(columns)]
    b = vector(SetUpMatrix([p],columns).list()[:-1])
    try:
        coeffs = M.solve_left(b)
    except:
        return False

    if not InputProd:
        cofactors = [(c*l.toNormal(),g.toNormal(),r.toNormal()) for (c,(l,g,r)) in zip(coeffs,cofactors) if c != 0]
    else:
        cofactors = [(NCMonomial(c*l.coeff,l.mon,intern=True),g,r) for (c,(l,g,r)) in zip(coeffs,cofactors) if c != 0]

    return cofactors
############################################################################
# Rewriting cofactors
############################################################################
def Rewrite(linear_comb,cofactors):
    certificate = []
    d = {f:c for (f,c) in cofactors}
    for (a,f,b) in linear_comb:
        try:
            cofactors_f = d[f]
            certificate += [(a*l,g,r*b) for (l,g,r) in cofactors_f]
        except:
            certificate.append((a,f,b))
    return certificate
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
def __symbolicPreprocessingRed__(h,G):
    T = {m for m in h.tail}
    T.add(h.lt)
    lt = [(g.lt.mon,g) for g in G]
    F = []
    columns = set()
    done = set()
    cofactors = []
    while len(T) > 0:
        t = T.pop()
        t_mon = t.mon
        columns.add(t)
        done.add(t_mon)
        for (mon,f) in lt:
            if mon in t_mon:
                factors = (t_mon).split(mon,1)
                l = factors[0] + '*'
                r = '*' + factors[1]
                g = f.lrmul(l,r)
                cofactors.append([NCMonomial(1,l,intern=True),f,NCMonomial(1,r,intern=True)])
                F.append(g)
                T.update({m for m in g.tail if m.mon not in done})
                break
    return [F,columns,cofactors]
############################################################################
def SetUpMatrix(F,columns):
    entries = {}
    l = len(columns)
    cols = {c.mon:i for (c,i) in zip(columns,range(l))}
    for i,f in enumerate(F):
        entries[(i,cols[f.lt.mon])] = f.lt.coeff
        entries[(i,i+l)] = 1
        entries.update({(i,cols[m.mon]):m.coeff for m in f.tail})

    return matrix(QQ,entries,sparse=True)
############################################################################
def Matrix2Poly(M,columns,F,cofactorsF):
    T = M[:,len(columns):] #transformation matrix
    M = M[:,:len(columns)]
    pos = M.nonzero_positions()
    #leave this line as is!
    FPlus = [[] for i in range(pos[-1][0]+1)]
    cofactors = [[] for i in range(pos[-1][0]+1)]
    for (i,j) in pos:
        FPlus[i].append(NCMonomial(M[i,j],columns[j].mon,intern=True))
    for (i,j) in T.nonzero_positions():
        if i < len(FPlus):
            f = copyTriple(cofactorsF[j])
            f[0].coeff *= T[i,j]
            cofactors[i].append(f)
        else:
            break
    lt_F = {f.lt.mon for f in F}
    FPlus = [(NCPoly(f[0],f[1:]),c) for (f,c) in zip(FPlus,cofactors) if f[0].mon not in lt_F]
 
    return FPlus
############################################################################
# Interreduction
############################################################################
def Interreduce(F,InputProd=False):
    if InputProd:
        G = F
    else:
        G = [NCPoly(f) for f in F]

    cofactors = [[(NCMonomial(1,'*',intern=True),g,NCMonomial(1,'*',intern=True))] for g in G]
    i = 0
    s = len(G)
    while i < s:
        if G[i] is None:
            i += 1
            continue
        (normal_form,linear_comb) = ReducedForm([g for (j,g) in enumerate(G) if j != i and g is not None],G[i],InputProd=True)

        #normal form is 0
        if type(normal_form) is int:
            G[i] = None
            cofactors[i] = None
        #normal form != 0
        elif not normal_form == G[i]:
            newPart = []
            for (a,f,b) in linear_comb:
                pos = G.index(f)
                newPart += [(a*l*(-1),h,r*b) for (l,h,r) in cofactors[pos]]
            cofactors[i] += newPart
            lc = normal_form.lt.coeff
            if lc != 1:
                normal_form.makeMonic()
                cofactors[i] = [(NCMonomial(a.coeff/lc,a.mon,intern=True),f,b) for (a,f,b) in cofactors[i]]
            G[i] = normal_form
            i = 0
        else:
            i += 1

    G = [g for g in G if g is not None]
    cofactors = [c for c in cofactors if c is not None]
    return (G,cofactors)
############################################################################
# Certify
############################################################################
def Certify(assumptionsInput,claims,Q,N=-1,MaxIter=10,MaxDeg=-1,MultiLex=False,Info=True,Criterion=True):
    
    assumptions = [NCPoly(f) for f in assumptionsInput]
    
    #checking compatibility
    for f in assumptions:
        if Q.QSignature(f) == []:
            print "The assumption " + str(f.toNormal()) + " is not compatible with the quiver."
            return False
    if type(claims) is list:
        for f in claims:
            if Q.QSignature(f) == []:
                print "The claim " + str(f.toNormal()) + " is not compatible with the quiver."
                return False
    else:
         if Q.QSignature(claims) == []:
            print "The claim " + str(claims.toNormal()) + " is not compatible with the quiver."
            return False
        
    #setting up the ring
    if Info:
        print "Using the following monomial ordering:"
    if MultiLex:
        vars_claims = set()
        if type(claims) is list:
            for f in claims:
                vars_claims.update(f.variables())
        else:
            vars_claims = claims.variables()
        vars_claims = [str(v) for v in vars_claims]
        knowns = [v for v in Q.vars() if v in vars_claims]
        unknowns = [v for v in Q.vars() if v not in knowns]
        SetUpRing(knowns,unknowns,Info=Info)
    else:
        SetUpRing(Q.vars(),Info=Info)

    assumptions = [NCPoly(f) for f in assumptionsInput]
    
    #keep track of leading coefficients
    lc = [f.lt.coeff for f in assumptions]
    for f in assumptions:
        f.makeMonic()

    #interreduce TODO
    (assumptionsRed,cofactorsRed) = Interreduce(assumptions,InputProd=True)
    if Info:
        print "\nInterreduced the input from " + str(len(assumptions)) + " polynomials to " + str(len(assumptionsRed)) + ".\n"

    #compute Groebner basis
    if Info:
        print "Computing a (partial) Groebner basis and reducing the claim...\n"
    cofactors = []
    i = 0
    linear_comb = false
    toIgnore = 0
    G = assumptionsRed

    while i < MaxIter and not linear_comb:
        toIgnoreOld = len(G)
        i+=1
        cofactors_temp = []
        if Info:
            print "Starting iteration " + str(i) + "...\n"
        G = F4(cofactors_temp,G,1,N=N,Ignore=toIgnore,Info=Info,MaxDeg=MaxDeg,Criterion=Criterion,OutputProd=True,IterCount=i-1)
        toIgnore = toIgnoreOld
        cofactors += cofactors_temp
        #try normal order
        if type(claims) is list:
            linear_comb = [IdealMembership(G,f,InputProd=True) for f in claims]
            if False in linear_comb:
                    linear_comb = False
        else:
            linear_comb = IdealMembership(G,claims,InputProd=True)

        #try random orders
        count = 0
        while count < 4 and not linear_comb:
            count += 1
            G_rand = sample(G,len(G))
            if type(claims) is list:
                linear_comb = [IdealMembership(G_rand,f,InputProd=True) for f in claims]
                if False in linear_comb:
                    linear_comb = False
            else:
                linear_comb = IdealMembership(G_rand,claims,InputProd=True)

    #negative outcome
    if not linear_comb:
        print "Not all claims could be reduced to 0."
        return False

    #positive outcome
    start = time()
    if Info:
        print "Rewriting the linear combination in terms of the assumptions..."
    certificate = __rewriteCertify__(linear_comb,cofactors)
    G = [f.toNormal() for f in G]
    #rewrite in terms of assumptions and not of the interreduced assumptions and take care of leading coefficients
    cofactorsRed = {g.toNormal():map(tripleToNormal,c) for (g,c) in zip(assumptionsRed,cofactorsRed)}
    lc = {f/c:c for (f,c) in zip(assumptionsInput,lc)}
    if type(claims) is list:
        for i in range(len(certificate)):
            cofactors_temp = []
            for (a,f,b) in certificate[i]:
                cofactors_temp += [(a*l,h,r*b) for (l,h,r) in cofactorsRed[f]]
            certificate[i] = [(a/lc[f],lc[f]*f,b) for (a,f,b) in cofactors_temp]
    else:
        cofactors_temp = []
        for (a,f,b) in certificate:
            cofactors_temp += [(a*l,h,r*b) for (l,h,r) in cofactorsRed[f]]
        certificate = [(a/lc[f],lc[f]*f,b) for (a,f,b) in cofactors_temp]
    end = time()
    if Info:
        print "Rewriting the linear combination took in total %.5f." % (end-start)
    if Info:
        print "\nDone! All claims were successfully reduced to 0."

    return certificate
############################################################################
def __rewriteCertify__(varsInput,cofactors):
    if type(varsInput[0]) is list:
        vars = flatten(varsInput)
    else:
        vars = varsInput

    #find all cofactors appearing in the linear combinations
    occurring = set()
    toAdd = set()
    for (l,f,r) in vars:
        toAdd.update({i for (i,(p,t)) in enumerate(cofactors) if f == p})
    while len(toAdd) > 0:
        occurring.update(toAdd)
        g = [cofactors[i] for i in toAdd]
        toAdd = set()
        for (h,t) in g:
            for (l,f,r) in t:
                toAdd.update({i for (i,(p,t)) in enumerate(cofactors) if f == p})
        toAdd.difference_update(occurring)

    if len(occurring) == 0:
        return map(tripleToNormal,vars)

    occurring = list(occurring)
    occurring.sort()

    #rewrite the appearing cofactors
    cofactors_done = {cofactors[occurring[0]][0].toNormal():map(tripleToNormal,cofactors[occurring[0]][1])}
    for i in occurring[1:]:
        cofactors_temp = []
        f = cofactors[i][0].toNormal()
        t = map(tripleToNormal,cofactors[i][1])
        for (a,g,b) in t:
            if g in cofactors_done:
                cofactors_temp += [(a*l,h,r*b) for (l,h,r) in cofactors_done[g]]
            else:
                cofactors_temp.append((a,g,b))
        cofactors_done[f] = cofactors_temp

    #rewrite the linear combinations
    if type(varsInput[0]) is list:
        l = len(varsInput)
        certificate = [[] for i in range(l)]
        for i in range(l):
            for (a,f,b) in map(tripleToNormal,varsInput[i]):
                if f in cofactors_done:
                    certificate[i] += [(a*l,g,r*b) for (l,g,r) in cofactors_done[f]]
                else:
                    certificate[i].append((a,f,b))
    else:
        certificate = []
        for (a,f,b) in map(tripleToNormal,varsInput):
            if f in cofactors_done:
                certificate += [(a*l,g,r*b) for (l,g,r) in cofactors_done[f]]
            else:
                certificate.append((a,f,b))

    return certificate
############################################################################
# Checking correctness
############################################################################
def MultiplyOut(cofactors):
    return sum(map(prod,cofactors))
############################################################################
def AllInI(cofactors,I):
    return all([g in I for [l,g,r] in cofactors])
############################################################################
def checkCorrect(cofactors,G,I):
    print "cofactors = G: "
    print [f for (f,c) in cofactors] == G[len(I):]

    print "linear combinations = G:"
    print map(MultiplyOut,[c for (f,c) in cofactors]) == G[len(I):]

    print "linear combination only of elements in I:"
    print all([AllInI(c,I) for (f,c) in cofactors])
############################################################################
def CheckCertificate(certificate,claim,I):
    print "certificate = claim:"
    print MultiplyOut(certificate) == claim

    print "certificate only of elements in I:"
    print all([g in I for (l,g,r) in certificate])
############################################################################
# additional stuff
############################################################################
def flatten(l):
    return [item for sublist in l for item in sublist]
############################################################################
def tripleToNormal(t):
    return (t[0].toNormal(),t[1].toNormal(),t[2].toNormal())
############################################################################
def copyTriple(t):
    return (t[0].copy(),t[1].copy(),t[2].copy())
############################################################################
#  Quiver
############################################################################
class Quiver:
    def __init__(self,triples):
        G = DiGraph(multiedges=True,loops=True)
        self.__vars__ = []
        for (l,s,t) in triples:
            G.add_edge(str(s),str(t),str(l))
            if l not in self.__vars__:
                self.__vars__.append(str(l))
        self.G = G
############################################################################
    def plot(self):
        G2 = self.G.plot(edge_labels=True,vertex_labels=False)
        G2.show()
############################################################################
    def vars(self):
        return self.__vars__
############################################################################
    def source(self,label):
        return [s for (s,t,l) in self.G.edge_iterator() if l == label]
############################################################################
    def target(self,label):
        return [t for (s,t,l) in self.G.edge_iterator() if l == label]
############################################################################
    def __QSignatureMon__(self,m):
        #if mon = 0 return V x V
        if m.coeff == 0:
            return {pair for pair in itertools.product(self.G.vertices(),repeat=2)}

        m = m.mon.strip('*').split('*')
        #if mon = const. return {(v,v) | v\in V}
        if m == ['']:
            return {(v,v) for v in self.G.vertex_iterator()}

        #usual case normal monomial
        begin = {(s,t) for (s,t,l) in self.G.edge_iterator() if l == m[-1]}
        for i in range(len(m)-2,-1,-1):
            end = {(s,t) for (s,t,l) in self.G.edge_iterator() if l == m[i]}
            comb = {(s1,t2) for ((s1,t1),(s2,t2)) in itertools.product(begin,end) if t1 == s2}
            if len(comb) == 0:
                return set()
            begin = comb
        return begin
############################################################################
    def QSignature(self,poly):
        if isinstance(poly,NCPoly):
            f = poly
        else:
            f = NCPoly(poly)
        monomials = f.tail + [f.lt]
        mon_signatures = [self.__QSignatureMon__(m) for m in monomials]
        return set.intersection(*mon_signatures)
 ############################################################################
    def compatibleQ(self,poly):
        return not len(self.QSignature(poly)) == 0
 ############################################################################
    def uniformlyCompatibleQ(self,poly):
        if not self.compatibleQ(poly):
            return False

        if isinstance(poly,NCPoly):
            f = poly
        else:
            f = NCPoly(poly)
        monomials = f.tail + [f.lt]
        mon_signatures = [self.__QSignatureMon__(m) for m in monomials]
        return all([sig == mon_signatures[0] for sig in mon_signatures])
############################################################################
# Stuff to execute
############################################################################
load("/Users/clemenshofstadler/Desktop/OperatorGB/Sage/Examples.py")
