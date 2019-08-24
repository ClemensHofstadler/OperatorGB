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
############################################################################
# polynomial datastructure
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
    
    def copy(self):
        return NCMonomial(self.coeff,copy(self.mon),intern=True)
    
    def __lt__(self,other):
        return SortedQ(self.mon,other.mon)

    def __repr__(self):
        return "(" + str(self.coeff) + ", " + self.mon + ")"
    
    def __eq__(self,other):
        return self.mon == other.mon and self.coeff == other.coeff
    
    def __hash__(self):
        return hash(self.mon)

    def __mul__(self,other):
        return NCMonomial(self.coeff*other.coeff,self.mon[:-1] + other.mon,intern = True)

    def toNormal(self):
        s = str(self.coeff)
        s += self.mon[:-1]

        F = FreeAlgebra(QQ,len(WordOrder),WordOrder)
        return F(s)


class NCPoly:
    def __init__(self,f,*args):
        if args:
            self.lt = f
            self.tail = args[0]
        else:
            d = f.monomial_coefficients()
            monomials = [NCMonomial(d[key],key) for key in d.keys()]
            monomials = sorted(monomials)
            self.lt = monomials[-1]
            self.tail = monomials[:-1]
    
    def copy(self):
        return NCPoly(self.lt.copy(),[mon.copy() for mon in self.tail])
    
    def __repr__(self):
        s = "(" + str(self.lt) + ", ["
        for m in self.tail:
            s += str(m) + ", "
            s = s[:-2] + "])"
        return s
    
    def __eq__(self,other):
        if len(self.tail) != len(other.tail):
            return False
        return self.lt == other.lt and all([m1 == m2 for (m1,m2) in zip(self.tail,other.tail)])
    
    def makeMonic(self):
        lc = self.lt.coeff
        if lc != 1:
            self.lt.coeff = 1
            for m in self.tail:
                m.coeff /= lc

    def lmul(self,m):
        f = self.copy()
        f.lt.mon = m[:-1] + f.lt.mon
        for mon in f.tail:
            mon.mon = m[:-1] + mon.mon
        return f

    def rmul(self,m):
        f = self.copy()
        f.lt.mon += m[1:]
        for mon in f.tail:
            mon.mon += m[1:]
        return f

    def toNormal(self):
        return self.lt.toNormal() + sum([_.toNormal() for _ in self.tail])
############################################################################
# set up the ring
############################################################################
def SetUpRing(vars,unknowns = None,Info=True):
    if unknowns != None:
        return __SetUpRing__(vars,unknowns)
    
    global WordOrder, SortedQ
    WordOrder = [str(v) for v in vars]
    SortedQ = DegLex
    if Info:
        s = ""
        for v in WordOrder:
            s += v + " < "
        print(s[:-3])

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
        print(s[:-3])
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
        amb = [a for a in amb if len(a.ABC.split('*'))-2 <= MaxDeg]

    return amb
############################################################################
def __GenerateAmbiguities__(oldPart,newPart,MaxDeg):
    
    amb = [generateOverlaps(v,w) + generateOverlaps(w,v) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [generateInclusions(v,w) + generateInclusions(w,v) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [generateOverlaps(v,w) for (v,w) in itertools.product(newPart,repeat=2)]
    amb += [generateInclusions(v,w) for (v,w) in itertools.product(newPart,repeat=2)]
    
    amb = flatten(amb)
    
    if MaxDeg > 0:
        amb = [a for a in amb if len(a.ABC.split('*'))-2 <= MaxDeg]

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
    amb.sort(key=lambda a:len(a.ABC.split('*')),reverse=True)
    indices = flatten([[a.i,a.j] for a in amb])
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
def CheckResolvability(G,oldlength=0,Criterion=True,MaxDeg=-1,Info=True,Sorted=False):
    words = [(f.lt.mon,i) for (i,f) in enumerate(G)]
    start = time()
    amb = GenerateAmbiguities(words[:oldlength],words[oldlength:],MaxDeg=MaxDeg)
    end = time()
    if Info:
        print(str(len(amb)) + " ambiguities in total (computation took %.5f)" % (end-start) )
    
    if Criterion:
        amb = DeleteRedundant(amb,Info=Info)
    if Sorted:
        print("Sorting the ambiguities is not implemented yet")

    spol = [SPoly(a,G[a.i],G[a.j]) for a in amb]
    spol = flatten(spol)
    if Info:
        print(str(len(spol)) + " critial polynoimals were generated.")

    return spol
############################################################################
def SPoly(amb,fi,fj):
    one = NCMonomial(1,'*',intern=True)
    mon_a = NCMonomial(1,amb.A,intern=True)
    mon_c = NCMonomial(1,amb.C,intern=True)
    if isinstance(amb,Inclusion):
        return [[fi,(one,fi,one)],
                [fj.lmul(amb.A).rmul(amb.C),(mon_a,fj,mon_c)]]
    else:
        return [[fi.rmul(amb.C),(one,fi,mon_c)],
                [fj.lmul(amb.A),(mon_a,fj,one)]]
############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
def Reduction(L,G):

    (L,cofactorsF) = map(list,zip(*L))
    (F,columns) = SymbolicPreprocessing(L,G,cofactorsF)
 
    columns = sorted(columns,reverse=True)
    M = SetUpMatrix(F,columns)
  
    M = M.rref()
    
    FPlus = Matrix2Poly(M,columns,F,cofactorsF)
    
    return FPlus
############################################################################
def SymbolicPreprocessing(L,G,cofactorsF):
    F = copy(L)
    T = {m for f in F for m in f.tail}
    lt = [(f.lt.mon,f) for f in G]
    columns = {f.lt for f in F}
    while len(T) > 0:
        t = T.pop()
        columns.add(t)
        for (mon,f) in lt:
            factors = (t.mon).split(mon,1)
            if len(factors) != 1:
                l = factors[0] + '*'
                r = '*' + factors[1]
                g = f.lmul(l).rmul(r)
                F.append(g)
                cofactorsF.append((NCMonomial(1,l,intern=True),f,NCMonomial(1,r,intern=True)))
                T.update({m for m in g.tail}-columns)
                break
    return [F,columns]
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
def rewriteCofactors(cofactors,I):

    cofactors[:] = [(f.toNormal(),map(tripleToNormal,c)) for (f,c) in cofactors]
    cofactors_done = {cofactors[0][0]:cofactors[0][1]}
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

def tripleToNormal(t):
    return (t[0].toNormal(),t[1].toNormal(),t[2].toNormal())

def copyTriple(t):
    return (t[0].copy(),t[1].copy(),t[2].copy())
############################################################################
# F4
############################################################################
def F4(cofactors,ideal,MaxIter=10,N=2,Ignore=0,MaxDeg=-1,Criterion=True,Info=True,Sorted=False,OutputProd=False):
    
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

    spol = CheckResolvability(G,Ignore,Criterion=Criterion,MaxDeg=MaxDeg,Info=Info,Sorted=Sorted)

    while count < MaxIter:
        start = time()
        while len(spol) > 0:
            n = min(N,len(spol))
            FPlus = Reduction(spol[:n],G)
            spol = spol[n:]
            F = [f for (f,c) in FPlus]
            G += F
            cofactors += FPlus
            print str(len(spol)).ljust(10)
            sys.stdout.write("\033[F")
   
        end = time()
        if Info:
            print "Reduction took: %.5f" % (end-start)
        if len(G) == oldlength:
            print "All S-polynomials could be reduced to 0."
            break
        
        count += 1
        if Info:
            print "Iteration " + str(count) + " finished. G has now " + str(len(G)) + " elements.\n"
        if count < MaxIter:
            spol = CheckResolvability(G,oldlength,Criterion=Criterion,MaxDeg=MaxDeg,Info=Info,Sorted=Sorted)
            oldlength = len(G)
            
    if OutputProd:
        return G
    else:
        start = time()
        if Info:
            print "Rewriting the cofactors has started."
        rewriteCofactors(cofactors,ideal)
        end = time()
        if Info:
            print "Rewriting the cofactors took in total %.5f\n" %(end-start)
        return [f.toNormal() for f in G]
            
############################################################################
# additional stuff
############################################################################
def flatten(l):
    return [item for sublist in l for item in sublist]
############################################################################
# Reduction
############################################################################
def SymbolicPreprocessingRed(h,G):
    T = {m for m in h.tail}
    T.add(h.lt)
    lt = [(g.lt.mon,g) for g in G]
    F = []
    columns = set()
    cofactors = []
    while len(T) > 0:
        t = T.pop()
        columns.add(t)
        for (mon,f) in lt:
            factors = (t.mon).split(mon,1)
            if len(factors) != 1:
                g = f.lmul(factors[0] + '*').rmul('*' + factors[1])
                cofactors.append([NCMonomial(1,factors[0] + '*',intern=True),f,NCMonomial(1,'*' + factors[1],intern=True)])
                F.append(g)
                T.update({m for m in g.tail}-columns)
                break
    return [F,columns,cofactors]

def ReducedForm(G_Input,f,InputProd=False):
    if not InputProd:
        G = map(NCPoly,G_Input)
    else:
        G = G_Input
    p = NCPoly(f)
    cofactors = []
    (F,columns,cofactors) = SymbolicPreprocessingRed(p,G)
    M = SetUpMatrix(F,columns)[:,:len(columns)]
    b = vector(SetUpMatrix([p],columns).list()[:-1])
    try:
        coeffs = M.solve_left(b)
    except:
        return False

    cofactors = [(c*l.toNormal(),g.toNormal(),r.toNormal()) for (c,(l,g,r)) in zip(coeffs,cofactors)]
    return cofactors

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
# Certify
############################################################################
def Certify(assumptions,claims,Q,N=2,MaxIter=10,MaxDeg=-1,MultiLex=False,Info=False,Criterion=True):
    
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
        knowns = [v for v in Q if v in vars_claims]
        unknowns = [v for v in Q if v not in knowns]
        SetUpRing(knowns,unknowns,Info=Info)
    else:
        SetUpRing(Q,Info=Info)

    G = [NCPoly(f) for f in assumptions]
    
    #keep track of leading coefficients
    lc = [f.lt.coeff for f in G]
    #lc = {f:c for (f,c) in zip(assumptions,lc)}

    #interreduce TODO

    #compute Groebner basis
    if Info:
        print "\nComputing a (partial) Groebner basis and reducing the claim...\n"
    cofactors = []
    i = 0
    linear_comb = false
    toIgnore = 0

    while i < MaxIter and not linear_comb:
        toIgnoreOld = len(G)
        i+=1
        cofactors_temp = []
        if Info:
            print "Starting iteration " + str(i) + "...\n"
        G = F4(cofactors_temp,G,1,N=N,Ignore=toIgnore,Info=Info,MaxDeg=MaxDeg,Criterion=Criterion,OutputProd=True)
        toIgnore = toIgnoreOld
        cofactors += cofactors_temp
        #try normal order
        if type(claims) is list:
            linear_comb = [ReducedForm(G,f,InputProd=True) for f in claims]
            if False in linear_comb:
                    linear_comb = False
        else:
            linear_comb = ReducedForm(G,claims,InputProd=True)

        #try random orders
        count = 0
        while count < 4 and not linear_comb:
            count += 1
            G_rand = sample(G,len(G))
            if type(claims) is list:
                linear_comb = [ReducedForm(G_rand,f,InputProd=True) for f in claims]
                if False in linear_comb:
                    linear_comb = False
            else:
                linear_comb = ReducedForm(G_rand,claims,InputProd=True)

    #negative outcome
    if not linear_comb:
        print "Not all claims could be reduced to 0."
        return False

    #positive outcome
    start = time()
    if Info:
        print "Rewriting the linear combination in terms of the assumptions..."
    rewriteCofactors(cofactors,G[:len(assumptions)])
    G = [f.toNormal() for f in G]
    #also take care of leading coefficients
    lc = {f/c:c for (f,c) in zip(assumptions,lc)}
    if type(claims) is list:
        certificate = [Rewrite(l,cofactors) for l in linear_comb]
        for i in range(len(certificate)):
            certificate[i] = [(a/lc[f],lc[f]*f,b) for (a,f,b) in certificate[i]]
    else:
        certificate = Rewrite(linear_comb,cofactors)
        certificate = [(a/lc[f],lc[f]*f,b) for (a,f,b) in certificate]
    end = time()
    if Info:
        print "Rewriting the linear combination took in total %.5f." % (end-start)
    if Info:
        print "\nDone! All claims were successfully reduced to 0."

    return certificate
############################################################################
# Checking correctness
############################################################################
def MultiplyOut(cofactors):
    return sum(map(prod,cofactors))

def AllInI(cofactors,I):
    return all([g in I for [l,g,r] in cofactors])

def checkCorrect(cofactors,G,I):
    print "cofactors = G: "
    print [f for (f,c) in cofactors] == G[len(I):]

    print "linear combinations = G:"
    print map(MultiplyOut,[c for (f,c) in cofactors]) == G[len(I):]

    print "linear combination only of elements in I:"
    print all([AllInI(c,I) for (f,c) in cofactors])

def CheckCertificate(certificate,claim,I):
    print "certificate = claim:"
    print MultiplyOut(certificate) == claim

    print "certificate only of elements in I:"
    print all([g in I for (l,g,r) in certificate])
############################################################################
# Stuff to execute
############################################################################
load("/Users/clemenshofstadler/Desktop/OperatorGB/Sage/Examples.py")
cofactors = []