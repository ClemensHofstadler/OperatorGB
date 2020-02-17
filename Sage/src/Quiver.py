# coding: utf-8
"""
Quiver
================

Class that allows to check (unifom) compatibility of polynomials with quivers.

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
from sage.all import DiGraph

from .NCMonomial import NCMonomial
from .NCPoly import NCPoly


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
        m = str(m).split('*')
        #if m = const. return {(v,v) | v\in V}
        if m == ['1']:
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
    def QSignature(self,f):
        monomials = f.monomials()
        if monomials:
            mon_signatures = [self.__QSignatureMon__(m) for m in monomials]
            return set.intersection(*mon_signatures)
        # f = 0
        else:
            return {pair for pair in itertools.product(self.G.vertices(),repeat=2)}
 ############################################################################
    def compatibleQ(self,f):
        return len(self.QSignature(f)) > 0
 ############################################################################
    def uniformlyCompatibleQ(self,f):
        if not self.compatibleQ(f):
            return False

        monomials = f.monomials()
        mon_signatures = [self.__QSignatureMon__(m) for m in monomials]
        return all([sig == mon_signatures[0] for sig in mon_signatures])
############################################################################
    def __repr__(self):
        vars = ""
        for v in self.vars():
            vars += v + ", "
        vars = vars[:-2]

        return "Labelled quiver with %s vertices in the labels {%s}." % (str(self.G.num_verts()),vars)
