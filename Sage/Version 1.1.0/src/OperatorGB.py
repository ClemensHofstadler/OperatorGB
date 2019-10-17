# coding: utf-8
"""
OperatorGB
================

Package to compute noncommutative Groebner bases in the free Algebra over QQ
together with tracing of cofactors. Additionally this package also allows to
check (unifom) compatibility of polynomials with quivers and allows to fully
automatically prove operator identities.

AUTHOR:

- Clemens Hofstadler (2019-07-08)

Version 1.1.0

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

from Ambiguities import flatten
from F4 import *
from NCMonomial import NCMonomial
from NCPoly import NCPoly
from Quiver import Quiver
from MonomialOrder import SetUpRing

############################################################################
# Certify
############################################################################
def Certify(assumptionsInput,claims,Q,N=-1,MaxIter=10,MaxDeg=-1,MultiLex=False,Info=True,Criterion=True):
    
    assumptions = [NCPoly(f) for f in assumptionsInput]
    
    #checking compatibility
    for f in assumptions:
        if not Q.uniformlyCompatibleQ(f):
            print "The assumption " + str(f.toNormal()) + " is not compatible with the quiver."
            return False
    if type(claims) is list:
        for f in claims:
            if not Q.compatibleQ(f):
                print "The claim " + str(f.toNormal()) + " is not compatible with the quiver."
                return False
    else:
        if not Q.compatibleQ(claims):
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

    #interreduction of the assumptions
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
        G = Groebner(cofactors_temp,G,1,N=N,Ignore=toIgnore,Info=Info,MaxDeg=MaxDeg,Criterion=Criterion,OutputProd=True,IterCount=i-1)
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

    print "certificate consists only of elements in I:"
    print all([g in I for (l,g,r) in certificate])
