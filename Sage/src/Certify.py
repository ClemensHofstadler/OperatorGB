# coding: utf-8
"""
Certify
================

Module to automatically prove operator identities.

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

from .MonomialOrder import SetUpRing
from .NCMonomial import NCMonomial
from .NCPoly import NCPoly
from .Ambiguities import flatten
from .NormalForms import __interreduce__,__idealMembership__
from .F4 import *
from .Quiver import Quiver

############################################################################
# Certify
############################################################################
def Certify(assumptionsInput,claimsInput,Q,MaxIter=10,MaxDeg=-1,MultiLex=False,Info=False,Criterion=True):
    
    # prepare input
    if type(claimsInput) is not list:
        claims = [claimsInput]
    else:
        claims = claimsInput
    
    # check compatibility
    __checkCompatibility__(assumptionsInput,claims,Q)
    
    # set up the ring
    __setUpRing__(claims,Q,Info,MultiLex)
    
    # prepare input further
    assumptions = [NCPoly(f) for f in assumptionsInput]
    claims = [NCPoly(f) for f in claims]
    for i,f in enumerate(assumptions):
        f.makeMonic(i)
    for f in claims:
        f.cofactors = False

    # interreduction of the assumptions
    G = __interreduce__(assumptions)
    if Info:
        print("\nInterreduced the input from " + str(len(assumptions)) + " polynomials to " + str(len(G)) + ".\n")

    # compute Groebner basis
    if Info:
        print("Computing a (partial) Groebner basis and reducing the claim...\n")
    i = 0
    toIgnore = 0
    keep_going = True
    cofactors = []

    while i < MaxIter and keep_going:
        toIgnoreOld = len(G)
        i += 1
        if Info:
            print("Starting iteration " + str(i) + "...\n")
        G = Groebner(cofactors,G,1,Ignore=toIgnore,Info=Info,MaxDeg=MaxDeg,Criterion=Criterion,IterCount=i-1)
        toIgnore = toIgnoreOld

        #try to reduce the claim using some (random) orders
        count = 0
        shuffle_flag = False
        while count < 4 and any(f.cofactors is False for f in claims):
            count += 1
            for f in claims:
                if f.cofactors is False:
                    f.cofactors = __idealMembership__(G,f,shuffle=shuffle_flag)
            shuffle_flag = True
        
        # if everything could be reduced => stop iteration
        keep_going = any(f.cofactors is False for f in claims)

    # negative outcome
    if any(f.cofactors == False for f in claims):
        print("Not all claims could be reduced to 0.")
        return False

    # positive outcome
    start = time()
    if Info:
        print("Rewriting the linear combination in terms of the assumptions...")
    certificate = __rewriteCertify__(claims,G,assumptionsInput)
    SimplifyCofactors(certificate)
    end = time()
    
    if Info:
        print("Rewriting the linear combination took in total %.5f." % (end-start))
        print("\nDone! All claims were successfully reduced to 0.")

    if type(claimsInput) is not list:
        return certificate[0]
    else:
        return certificate
############################################################################
def __checkCompatibility__(assumptions,claims,Q):
    for f in assumptions:
        if not Q.uniformlyCompatibleQ(f):
            raise ValueError("The assumption " + str(f.toNormal()) + " is not compatible with the quiver.")
    for f in claims:
        if not Q.compatibleQ(f):
                raise ValueError("The claim " + str(f) + " is not compatible with the quiver.")
############################################################################
def __setUpRing__(claims,Q,Info,MultiLex):
    if Info:
        print("Using the following monomial ordering:")
    if MultiLex:
        vars_claims = set()
        for f in claims:
            vars_claims.update(f.variables())

        vars_claims = [str(v) for v in vars_claims]
        knowns = [v for v in Q.vars() if v in vars_claims]
        unknowns = [v for v in Q.vars() if v not in knowns]
        SetUpRing(knowns,unknowns,Info=Info)
    else:
        SetUpRing(Q.vars(),Info=Info)
############################################################################
def __rewriteCertify__(claims,G,F):
    # find all elements of G that occur in the linear combinations
    N = len(F)
    P = F[0].parent()
    occurring = set()
    done = set()
    for f in claims:
        occurring.update({i for (a,i,b) in f.cofactors})
        f.cofactors = [__tripToNorm__((a,j,b),P) for (a,j,b) in f.cofactors]
    while len(occurring) > 0:
        i = occurring.pop()
        done.add(i)
        occurring.update({j for (a,j,b) in G[i].cofactors if j not in done})
    
    occurring = done.union(set(range(N)))
    # rewrite the cofactors of the occuring elements in G
    G_occurring = {i:g for i,g in enumerate(G) if i in occurring}
    occurring = list(occurring)
    occurring.sort()
    cofactors = {i:[(a.toNormal(P),j,b.toNormal(P)) for (a,j,b) in G_occurring[i].cofactors] for i in occurring}
    for i in range(N):
        cofactors[i] = [(a,F[j],b) for (a,j,b) in cofactors[i]]
    for j in occurring[N:]:
        cofactors_temp = []
        for (a,i,b) in cofactors[j]:
            cofactors_temp += [(a*l,f,r*b) for (l,f,r) in cofactors[i]]
        cofactors[j] = cofactors_temp

    # rewrite the cofactors of the claims
    certificate = [[] for i in range(len(claims))]
    for idx,f in enumerate(claims):
        for (a,i,b) in f.cofactors:
            certificate[idx] += [(a*l,g,r*b) for (l,g,r) in cofactors[i]]
 
    return certificate
############################################################################
def __tripToNorm__(t,P):
    return (t[0].toNormal(P),t[1],t[2].toNormal(P))
