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




