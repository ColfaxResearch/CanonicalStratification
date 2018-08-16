""" Script for bencmarking computational
performance of stratification algorith for
the 2-sphere and 3-ball.

Distributed under MIT License

Authors: Ryo Asai (ryo@colfax-intl.com)
         Jay Shah (jshah3@nd.edu)
"""
import sys
from CanonicalStratification.CanonicalStratification import *
from CanonicalStratification.Utils import *

# Required for connected component analysis
sys.setrecursionlimit(1000000)


# Available method for constructing the simplicial complex
complex_constructors={"ball"    : ballDelaunay,
                      "2sphere" : sphereDelaunay}

if(len(sys.argv) < 3):
  print("Usage: "+sys.argv[0]+" complex_type vertex_count\n  complex_type is one of: "+' '.join(complex_constructors.keys()))
  exit(1)
complex_type = sys.argv[1];
vertex_count = sys.argv[2];

if complex_type in complex_constructors: 
  simplexes = complex_constructors[complex_type](int(vertex_count));
  sc = SimplicialComplex.fromList(simplexes)
  print("%d simplexes in the complex" % sc.getNumberOfSimplexes())
  for i in range(10):
    sm = stratify(sc)
    del sm
else:
  print("Unknown comlplex type "+complex_type)
  print(" available types are: "+str(complex_constructors.keys()))
