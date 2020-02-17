# coding: utf-8
"""
Quiver
================

Class that allows to check (unifom) compatibility of polynomials with quivers.

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
import itertools
from sage.all import DiGraph

from NCMonomial import NCMonomial
from NCPoly import NCPoly


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
    def __repr__(self):
        vars = ""
        for v in self.vars():
            vars += v + ", "
        vars = vars[:-2]

        return "Labelled quiver with %s vertices in the labels {%s}." % (str(self.G.num_verts()),vars)
