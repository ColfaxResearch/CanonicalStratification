""" Script for bencmarking computational
performance of stratification algorith for
the 2-sphere and 3-ball.

(c) Colfax Research
Distributed under MIT License

Authors: Ryo Asai (ryo@colfax-intl.com)
         Jay Shah (jshah3@nd.edu)
"""
import sys
from CanonicalStratification.CanonicalStratification import *

from scipy.spatial import Delaunay
import numpy as np

def ballDelaunay(num_vertexes):
  ''' Produces a triangulation of 3-ball. '''
  sph_coord = np.random.random_sample((3,num_vertexes))*np.array([1.,2.*np.pi,np.pi]).reshape(3,1)
  ecl_coord = np.empty_like(sph_coord)
  ecl_coord[0,:]=sph_coord[0,:]*np.cos(sph_coord[1,:])*np.sin(sph_coord[2,:])
  ecl_coord[1,:]=sph_coord[0,:]*np.sin(sph_coord[1,:])*np.sin(sph_coord[2,:])
  ecl_coord[2,:]=sph_coord[0,:]*np.cos(sph_coord[2,:])
  return Delaunay(ecl_coord.T).simplices

def sphereDelaunay(num_vertexes):
  ''' Produces a triangulation of 2-sphere. '''
  ang_coord = np.random.random_sample((2,num_vertexes))*np.array([2.*np.pi,np.pi]).reshape(2,1)
  ecl_coord = np.empty((3,num_vertexes), dtype=np.float)
  ecl_coord[0,:]=np.cos(ang_coord[0,:])*np.sin(ang_coord[1,:])
  ecl_coord[1,:]=np.sin(ang_coord[0,:])*np.sin(ang_coord[1,:])
  ecl_coord[2,:]=np.cos(ang_coord[1,:])
  return Delaunay(ecl_coord.T).simplices




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
