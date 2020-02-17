# coding: utf-8
"""
MonomialOrder
================

Module to set different monomial orderings for Greobner basis computations
in the free algebra.

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
import GlobalData

############################################################################
# set up the ring
############################################################################
def SetUpRing(vars,unknowns = None,Info=True):
    if unknowns != None:
        return __SetUpRing__(vars,unknowns,Info=Info)
    
    GlobalData.WordOrder = [str(v) for v in vars]
    GlobalData.MonomialOrder = "DegLex"
    if Info:
        s = ""
        for v in GlobalData.WordOrder:
            s += v + " < "
        print s[:-3]
############################################################################
def __SetUpRing__(knownsInput,unknownsInput,Info=True):
   
    knowns = [str(v) for v in knownsInput]
    unknowns = [str(v) for v in unknownsInput]
    GlobalData.WordOrder = knowns + unknowns
    GlobalData.Knowns = knowns
    GlobalData.Unknowns = unknowns
    GlobalData.MonomialOrder = "MultiLex"
    if Info:
        s = ""
        for v in knowns:
            s += v + " < "
        s = s[:-1]
        s += "< "
        for v in unknowns:
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
            return GlobalData.WordOrder.index(a_vars[i]) < GlobalData.WordOrder.index(b_vars[i])

    else:
        return la < lb
############################################################################
def MultiLex(a,b):
    knowns = GlobalData.Knowns
    unknowns = GlobalData.Unknowns
    
    a_vars = a.split('*')
    b_vars = b.split('*')
    V1a = len([v for v in a_vars if v in unknowns])
    V2a = len([v for v in a_vars if v in knowns])
    V1b = len([v for v in b_vars if v in unknowns])
    V2b = len([v for v in b_vars if v in knowns])

    if V1a < V1b or (V1a == V1b and V2a < V2b):
        return True
    if V1a > V1b or (V1a == V1b and V2a > V2b):
        return False

    return DegLex(a,b)
############################################################################
# SortedQ
############################################################################
def SortedQ(a,b):
    if GlobalData.MonomialOrder == "DegLex":
        return DegLex(a,b)
    elif GlobalData.MonomialOrder == "MultiLex":
        return MultiLex(a,b)
    else:
        return None
