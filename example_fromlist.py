""" Example script for loading simplices
from a file.
 
Distributed under MIT License

Authors: Ryo Asai (ryo@colfax-intl.com)
         Jay Shah (jshah3@nd.edu)
"""
import sys
from CanonicalStratification import *

simplices=[[0,1,3],[0,2,3],[1,3,5],[2,3,4],[2,4,6],[3,4,5],[4,5,7],[4,6,7],[3,4,8]]

sc = SimplicialComplex.fromList(simplices)

sm = stratify(sc)
print("Stratification found %d strata" % (sm.strata_count,))

strata_members = [set() for _ in range(sm.strata_count)]
for d_simplices in sc.simplices:
  for simplex in d_simplices.values():
    strata_members[sm.getStratum(simplex).index].add(simplex)

for index in range(len(strata_members)):
  print("Stratum %d: %s " % (index,strata_members[index]))

