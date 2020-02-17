# coding: utf-8
"""
MonomialOrder
================

Module to set different monomial orderings for Greobner basis computations
in the free algebra.

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

from . import GlobalData

############################################################################
# set up the ring
############################################################################
def SetUpRing(*argv,Info=True):
    
    if len(argv) > 1:
        return __SetUpMultiLex__(*argv,Info=Info)
    
    GlobalData.WordOrder = [str(v) for v in argv[0]]
    GlobalData.WordOrderDict = {v:i for i,v in enumerate(GlobalData.WordOrder)}
    GlobalData.MonomialOrder = "DegLex"

    if Info:
        s = ""
        for v in GlobalData.WordOrder:
            s += v + " < "
        print(s[:-3])
############################################################################
def __SetUpMultiLex__(*argv,Info=True):
   
    GlobalData.MonomialOrder = "MultiLex"
    varSets = []
    for arg in argv:
        varSets.append([str(v) for v in arg])
    GlobalData.WordOrder = flatten(varSets)
    GlobalData.WordOrderDict = {v:i for i,v in enumerate(GlobalData.WordOrder)}
    
    if Info:
        s = ""
        for var in varSets:
            if len(var) > 0:
                for v in var[:-1]:
                    s += v + " < "
                s += var[-1] + " << "
        print(s[:-3])

    varSets.reverse()
    for i,varSet in enumerate(varSets):
        varSets[i] = set(varSet)
    GlobalData.VarSets = varSets
############################################################################
# monomial order
############################################################################
def DegLex(a,b):
    a_vars = a.split('*')[1:-1]
    b_vars = b.split('*')[1:-1]
    
    la = len(a_vars)
    lb = len(b_vars)

    if la == lb:
        x,y = next(((x,y) for (x,y) in zip(a_vars, b_vars) if x!=y),(None,None))
        if x:
            return GlobalData.WordOrderDict[x] < GlobalData.WordOrderDict[y]
        else:
            return True
    else:
        return la < lb
############################################################################
def MultiLex(a,b):
    a_vars = a.split('*')[1:-1]
    b_vars = b.split('*')[1:-1]

    for varSet in GlobalData.VarSets:
        Va = sum(v in varSet for v in a_vars)
        Vb = sum(v in varSet for v in b_vars)
        if Va < Vb:
            return True
        if Va > Vb:
            return False

    return __degLex__(a_vars,b_vars)
############################################################################
def __degLex__(a_vars,b_vars):
    x,y = next(((x,y) for (x,y) in zip(a_vars, b_vars) if x!=y),(None,None))
    if x:
        return GlobalData.WordOrderDict[x] < GlobalData.WordOrderDict[y]
    else:
        return True
############################################################################
# SortedQ
############################################################################
def SortedQ(a,b):
    if GlobalData.MonomialOrder == "DegLex":
        return DegLex(a,b)
    elif GlobalData.MonomialOrder == "MultiLex":
        return MultiLex(a,b)
    else:
        raise NotImplementedError
############################################################################
# additional stuff
############################################################################
def flatten(l):
    return [item for sublist in l for item in sublist]
