"""
OperatorGB
================

Package to compute noncommutative Groebner bases in the free Algebra over QQ
together with tracing of cofactors. Additionally this package also allows to
check (unifom) compatibility of polynomials with quivers and allows to fully
automatically prove operator identities.

AUTHOR:

- Clemens Hofstadler (2020-02-16)

Version 1.2.0

"""

#############################################################################
#  Copyright (C) 2020 Clemens Hofstadler (clemens.hofstadler@jku.at).       #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .src.MonomialOrder import SetUpRing,SortedQ
from .src.F4 import Groebner
from .src.NormalForms import ReducedForm,Rewrite,MultiplyOut,CheckCofactors,CheckCertificate
from .src.Certify import Certify
from .src.Quiver import Quiver

############################################################################
# Info
############################################################################
print("Package OperarorGB version 1.2.0")
print("Copyright 2020, Institute for Algebra, JKU")
print("by Clemens Hofstadler, clemens.hofstadler@jku.at")
