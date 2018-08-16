""" Example script for loading simplices
from a file.
 
Distributed under MIT License

Authors: Ryo Asai (ryo@colfax-intl.com)
         Jay Shah (jshah3@nd.edu)
"""
import sys
from CanonicalStratification import *

sys.setrecursionlimit(1000000)

sc = SimplicialComplex.fromFile("genus-50-surface.txt")

sm = stratify(sc)
print("Stratification found %d strata" % (sm.strata_count,))
test_simplex = [0,522,1408]
stratum =sm.getStratum(test_simplex)
print("Simplex %s is in stratum %d (dimension %d)" % (test_simplex, stratum.index, stratum.top_dim))
