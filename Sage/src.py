from __future__ import absolute_import

import itertools
from time import time
from random import choice
from sage.all import *
from sage.matrix.matrix_rational_sparse import Matrix_rational_sparse
from collections import defaultdict, deque

WordOrder = []
VarSets = []
MonomialOrder = None

"""
TODO when splitting into several classes:
 - remove "global" from SetUpRing functions

"""

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
                        mon += v + '*'
  
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
    def one():
        return NCMonomial(1,'*',intern=True)
############################################################################
    def zero():
        return NCMonomial(0,'1')
############################################################################
    def toNormal(self,Parent=None):
        s = str(self.coeff)
        s += self.mon[:-1]
        if not Parent:
            Parent = FreeAlgebra(QQ,len(WordOrder),WordOrder)
        return Parent(s)
########################################################################################################################################################
########################################################################################################################################################
# NCPoly
########################################################################################################################################################
########################################################################################################################################################
class NCPoly:
    def __init__(self,f,*args):
        if args:
            self.lt = f
            self.tail = args[0]
            if len(args) > 1:
                self.cofactors = args[1]
        else:
            d = f.monomial_coefficients()
            if d:
                monomials = [NCMonomial(d[key],key) for key in d.keys()]
                monomials = sorted(monomials)
            else:
                monomials = [NCMonomial(0,'*',intern=True)]
            self.lt = monomials[-1]
            self.tail = monomials[:-1]
            
            self.cofactors = None
        self.redundant = False
############################################################################
    def copy(self):
        return NCPoly(self.lt.copy(),[mon.copy() for mon in self.tail])
############################################################################
    def zero():
        return NCPoly(NCMonomial.zero(),[],[])
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
        if other is None:
            return False
        if len(self.tail) != len(other.tail):
            return False
        return self.lt == other.lt and all([m1 == m2 for (m1,m2) in zip(self.tail,other.tail)])
############################################################################
    def __hash__(self):
        return hash((self.lt,) + tuple(self.tail))
############################################################################
    def makeMonic(self,i = None):
        """
        Really update self
        """
        lc = self.lt.coeff
        if lc != 1:
            self.lt.coeff = 1
            for m in self.tail:
                m.coeff /= lc
        if i is not None:
            self.cofactors = [(NCMonomial(1/lc,"*",intern=True),i,NCMonomial.one())]
############################################################################
    def lmul(self,m):
        """
        Only update copy
        """
        f = self.copy()
        f.lt.mon = m[:-1] + f.lt.mon
        for mon in f.tail:
            mon.mon = m[:-1] + mon.mon
   
        return f
############################################################################
    def rmul(self,m):
        """
        Only update copy
        """
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
    def toNormal(self,Parent=None):
        return self.lt.toNormal(Parent) + sum([_.toNormal(Parent) for _ in self.tail])
############################################################################
# set up the ring
############################################################################
def SetUpRing(*argv,Info=True):
    global WordOrder,MonomialOrder,WordOrderDict
    
    if len(argv) > 1:
        return __SetUpMultiLex__(*argv,Info=Info)
    
    WordOrder = [str(v) for v in argv[0]]
    WordOrderDict = {v:i for i,v in enumerate(WordOrder)}
    MonomialOrder = "DegLex"

    if Info:
        s = ""
        for v in WordOrder:
            s += v + " < "
        print(s[:-3])
############################################################################
def __SetUpMultiLex__(*argv,Info=True):
    global WordOrder,MonomialOrder,VarSets,WordOrderDict
    
    MonomialOrder = "MultiLex"
    VarSets = []
    for arg in argv:
        VarSets.append([str(v) for v in arg])
    WordOrder = flatten(VarSets)
    WordOrderDict = {v:i for i,v in enumerate(WordOrder)}
    
    if Info:
        s = ""
        for var in VarSets:
            if len(var) > 0:
                for v in var[:-1]:
                    s += v + " < "
                s += var[-1] + " << "
        print(s[:-3])

    VarSets.reverse()
    for i,varSet in enumerate(VarSets):
        VarSets[i] = set(varSet)
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
            return WordOrderDict[x] < WordOrderDict[y]
        else:
            return True
    else:
        return la < lb
############################################################################
def MultiLex(a,b):
    a_vars = a.split('*')[1:-1]
    b_vars = b.split('*')[1:-1]

    for varSet in VarSets:
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
        return WordOrderDict[x] < WordOrderDict[y]
    else:
        return True
############################################################################
# SortedQ
############################################################################
def SortedQ(a,b):
    if MonomialOrder == "DegLex":
        return DegLex(a,b)
    elif MonomialOrder == "MultiLex":
        return MultiLex(a,b)
    else:
        raise NotImplementedError

############################################################################
# ambiguities
############################################################################
class Overlap:
    """
    fiC - Afj
    """
    def __init__(self,ABC,A,C,i,j):
        self.ABC = ABC
        self.A = A
        self.C = C
        self.i = i
        self.j = j
        self.deg = len(ABC.split('*'))
        self.min = min(i,j)

    def __repr__(self):
        return "Overlap(" + str(self.ABC) + ", " + str(self.A) + ", " + str(self.C) + ", (" + str(self.i) + ", " + str(self.j) + "))"
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
        self.min = min(i,j)

    def __repr__(self):
        return "Inclusion(" + str(self.ABC) + ", " + str(self.A) + ", " + str(self.C) + ", (" + str(self.i) + ", " + str(self.j) + "))"
############################################################################
def GenerateAmbiguities(words, newPart = None, MaxDeg = -1):
    if newPart != None:
        return __generateAmbiguities__(words,newPart,MaxDeg)
    
    amb = [generateOverlaps(v,w) for (v,w) in itertools.product(words,repeat=2)]
    amb += [list(generateInclusions(v,w)) for (v,w) in itertools.product(words,repeat=2)]

    amb = flatten(amb)
   
    if MaxDeg > 0:
        amb = [a for a in amb if a.deg-2 <= MaxDeg]

    return amb
############################################################################
def __generateAmbiguities__(oldPart,newPart,MaxDeg):
    
    amb = [generateOverlaps(v,w) + generateOverlaps(w,v) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [list(generateInclusions(v,w)) + list(generateInclusions(w,v)) for (v,w) in itertools.product(oldPart,newPart)]
    amb += [generateOverlaps(v,w) for (v,w) in itertools.product(newPart,repeat=2)]
    amb += [list(generateInclusions(v,w)) for (v,w) in itertools.product(newPart,repeat=2)]
    
    amb = flatten(amb)
    
    if MaxDeg > 0:
        amb = [a for a in amb if a.deg-2 <= MaxDeg]

    return amb
############################################################################
def generateOverlaps(a,b):
    v = a[0]
    w = b[0]
    return [Overlap(v + w[k:],v[:-k+1],w[k-1:],a[1],b[1]) for k in range(2,min(len(v),len(w))-1) if v.endswith(w[:k])]
############################################################################
def generateInclusions(a,b):
    v = a[0]
    w = b[0]
    k = 0
    if len(w) <= len(v) and a[1] != b[1]:
        while True:
            k = v.find(w, k)+1
            if k > 0:
                yield Inclusion(v,v[:k],v[k+len(w)-2:],a[1],b[1])
            else:
                break
############################################################################
def DeleteRedundant(amb,lt,Info=True):
    start = time()
    overlaps = [a for a in amb if type(a) == Overlap]
    incls = [a for a in amb if type(a) == Inclusion]
    
    for V,k in lt:
        overlaps = [a for a in overlaps if not (k < a.min and V in a.ABC)]
        incls = [a for a in incls if not (k < a.j and _incl_test_(a,V))]
    
    result = overlaps + incls
    if Info:
        print("Removed " + str(len(amb) - len(result)) + " ambiguities in %.5f" % (time()-start))
    return result

def _incl_test_(a,V):
    """
    V|A or V|B or V|C or B|V|ABC but V != ABC
    """
    if V in a.A or V in a.C:
        return True
    B = a.ABC[len(a.A)-1 : -len(a.C)+1]
    if V in B or (B in V and V in a.ABC[1:-1]):
        return True
    else:
        return False
############################################################################
# additional stuff
############################################################################
def flatten(l):
    return [item for sublist in l for item in sublist]
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
def SetUpMatrix(F,columns):
    entries = {}
    l = len(columns)
    cols = {c.mon:i for i,c in enumerate(columns)}
    for i,f in enumerate(F):
        entries[(i,cols[f.lt.mon])] = f.lt.coeff
        entries[(i,i+l)] = 1
        entries.update({(i,cols[m.mon]):m.coeff for m in f.tail})

    return matrix(QQ,entries,sparse=True)
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
def __tripToNorm__(t,P):
    return (t[0].toNormal(P),t[1],t[2].toNormal(P))
############################################################################
def copyTriple(t):
    return (t[0].copy(),t[1],t[2].copy())
############################################################################
def multiplyTriple(triple,c):
    t = copyTriple(triple)
    t[0].coeff *= c
    return t
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
# Normal form computation
############################################################################
def __interreduce__(G):
    i = 0
    s = len(G)
    while i < s:
        if G[i] is None:
            i += 1
            continue
    
        normal_form = __reducedForm__(G,i)
        # normal form = 0
        if type(normal_form) is int:
            G[i] = None
        # normal form != 0
        elif normal_form:
            G[i] = normal_form
            i = 0
        # normal form = g_i
        else:
            i += 1

    return [g for g in G if g]
############################################################################
def ReducedForm(cofactorsF,G_input,f_input):
    G = [NCPoly(g) for g in G_input]
    for i,g in enumerate(G):
        g.makeMonic(i)
    f = NCPoly(f_input)
    f.cofactors = [(NCMonomial.one(),len(G),NCMonomial.one())]

    normal_form = __reducedForm__(G+[f],len(G),intern=False)

    cofactorsF[:] = [tripleToNormal((a.__mul__(-1),G[i],b)) for (a,i,b) in normal_form.cofactors[1:]]
    return normal_form.toNormal()

def Rewrite(cofactorsF, cofactorsG):
    G = list(map(MultiplyOut,cofactorsG))
    result = []
    for (a,g,b) in cofactorsF:
        i = G.index(g)
        result += [(a*l,f,r*b) for (l,f,r) in cofactorsG[i]]
    return result

############################################################################
def __reducedForm__(G,i,intern=True):
    f = G[i]
    F = [f]
    columns,cofactors = __symbolicPreprocessingRed__(F,G,i)
    cofactors = [(NCMonomial.one(),i,NCMonomial.one())] + cofactors
    columns = [NCMonomial(1,c,intern=True) for c in columns]
    columns.sort(reverse=True)

    M = SetUpMatrix(F,columns).rref()
    T = M[:,len(columns):] #transformation matrix
    RREF = M[:,:len(columns)]

    row_idx = T.nonzero_positions_in_column(0)[-1]

    # reduction to zero
    if len(RREF.nonzero_positions_in_row(row_idx)) == 0:
        if intern:
            return 0
        else:
            normal_form = NCPoly.zero()
    else:
        monomials = []
        for j in RREF.nonzero_positions_in_row(row_idx):
            monomials.append(NCMonomial(RREF[row_idx,j],columns[j].mon,intern=True))
        normal_form = NCPoly(monomials[0],monomials[1:],[])

    # no reduction
    if intern and normal_form == f:
        return None
            
    # some reduction - take care of cofactors
    coeff = T[row_idx,0]
    for j in T.nonzero_positions_in_row(row_idx):
        (a,idx,b) = cofactors[j]
        coeff = T[row_idx,j]
        normal_form.cofactors += [(a*l*coeff,k,r*b) for (l,k,r) in G[idx].cofactors]
    return normal_form

############################################################################
def __symbolicPreprocessingRed__(F,G,idx,shuffle=False):
    T = {m.mon for m in F[0].tail}
    T.add(F[0].lt.mon)
    lt = [(i,g.lt.mon) for i,g in enumerate(G) if i != idx and g is not None]
    done = set()
    cofactors = []
    while len(T) > 0:
        t = T.pop()
        done.add(t)
        if shuffle:
            choices = [(i,m) for i,m in lt if m in t]
            i,m = choice(choices) if choices else (None,None)
        else:
            i,m = next(((i,m) for i,m in lt if m in t), (None,None))
        if m:
            g = G[i]
            factors = t.split(m,1)
            a = factors[0] + '*'
            b = '*' + factors[1]
            agb = g.lrmul(a,b)
            cofactors.append([NCMonomial(1,a,intern=True),i,NCMonomial(1,b,intern=True)])
            F.append(agb)
            T.update({m.mon for m in agb.tail if m.mon not in done})
    return done, cofactors
    ############################################################################
def __idealMembership__(G,f,shuffle=False):
    F = [f]
    columns,cofactors = __symbolicPreprocessingRed__(F,G,-1,shuffle=shuffle)
    columns = [NCMonomial(1,c,intern=True) for c in columns]
    columns.sort(reverse=True)
    M = SetUpMatrix(F[1:],columns)[:,:len(columns)]
    b = vector(SetUpMatrix(F[:1],columns).list()[:-1])
    try:
        coeffs = M.solve_left(b)
    except:
        return False

    cofactors = [multiplyTriple(t,c) for c,t in zip(coeffs,cofactors) if c != 0]

    return cofactors
############################################################################
# Interreduction
############################################################################

############################################################################
# Checking correctness
############################################################################
def MultiplyOut(cofactors):
    return sum(map(prod,cofactors))
############################################################################
def AllInI(cofactors,I):
    return all([g in I for (l,g,r) in cofactors])
############################################################################
def CheckCofactors(cofactors,G,I):
    print("cofactors = G:")
    print(all(MultiplyOut(c) == G[i] for i,c in enumerate(cofactors)))

    print("linear combination only of elements in I:")
    print(all(AllInI(c,I) for c in cofactors))
############################################################################
def CheckCertificate(certificate,claim,I):
    print("certificate = claim:")
    print(MultiplyOut(certificate) == claim)

    print("certificate consists only of elements in I:")
    print(all(g in I for (l,g,r) in certificate))
############################################################################
def CheckInterreduction(I,G):
    P = G[0].parent()
    cofactors = []
    for f in I:
        cofactors.append([(a.toNormal(P),G[j],b.toNormal(P)) for (a,j,b) in f.cofactors])
    print(all(MultiplyOut(c) == f.toNormal(P) for c,f in zip(cofactors,I)))
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
############################################################################
# Info
############################################################################
print("Package OperarorGB version 1.2.0")
print("Copyright 2019, Institute for Algebra, JKU")
print("by Clemens Hofstadler, clemens.hofstadler@jku.at")
