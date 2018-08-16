""" Class for containing the simplicial complex

Distributed under MIT License

Authors: Ryo Asai (ryo@colfax-intl.com)
         Jay Shah (jshah3@nd.edu)
"""
import itertools

class SimplicialComplex:
  ''' Class containing the Simplicial Complex
  
  This class is used as the input to the 
  canonical stratification. Currently, we
  support loading from an iterable of lists,
  (e.g. list of lists) and a from file.

  - SimplicialComplex(iterable_of_lists):
  - SimplicialComplex.fromList(simplex_list):
  - SimplcicalComplex.fromFile(filename):

  See the respective constructors for examples.

  The complex is a multi-partite graph with
  simplices as graph nodes. Each simplex node
  only has edges to its immediate coface or
  immediate face. Edge itself is not represented
  in this graph.
  
  To allow for easier access, the smplices are
  stored in a list of dictionaries. Each dict
  conatins all simplices of the same dimesnion.
  The key for the dict is the sorted of vertices
  of the simplex, which is unique. 
  '''


  class SimplexNode:
    ''' Class for a node in the graph. Represents
    a simplex.

    Should not be manually created, but instead
    through addSimplex()  
    '''
    def __init__(self, vertices, hash_val):
      # vertices making up the simplex
      self.vertices=vertices
      self.dim=len(vertices)-1
      self.hash_val = hash_val

      # "edges" to immediate coface/face
      self.cofaces=[]
      self.faces=[]

    def __str__(self):
      return str(self.vertices)
  
    def __repr__(self):
      return str(self)

  def addSimplex(self, vertices):
    ''' Function for adding a simplex to the 
    complex

    Works like a set in that it only adds the 
    simplex if it does not 

    return simplex,new

    - simplex - SimplexNode object corresponding 
                to the vertices
    - new     - Boolean; True if it is new. False 
                otherwise.
    '''

    # Check that the dict exists for this deimension
    dim = len(vertices)-1
    if len(self.simplices) -1 < dim:
      # Starting a new dict for this dimension
      while len(self.simplices) -1 < dim:
        self.simplices.append({})

    # Getting the dictionary key
    hash_val = self.getSimplexHash(vertices)

    # Check if the simplex already exists
    new_face =hash_val not in self.simplices[dim]
    if new_face:
      new_s = self.SimplexNode(vertices, hash_val)
      self.simplices[dim][hash_val]=new_s
    return self.simplices[dim][hash_val], new_face

  @staticmethod
  def getSimplexHash(vertices):
    ''' Returns the "hash" from a vertex set.
    
    The hash value is just the sorted of vertices 
    of the simplex, which is unique. 
    This is generally used as the key for various
    dict objects.
    '''
    return str(sorted([int(i) for i in vertices]))

  def __init__(self, iterable_simplices):
    ''' Constructs the complex from iterable of lists

    Each item in the iterable should contain an integer
    list with the indexes of the vertices in the simplex.
    Example:
       iterable_simplices = [[0,1,2],[0,3]]
    
    By definition, each subset of the simplex must also
    be in the complex, so we create all the implied 
    simplices. For example, [[0,1,2],[0,3]] will create:
    [[0],[1],[2],[3],[0,1],[1,2],[0,2],[0,3],[1,2,3]]
    '''

    self.simplices=[{}]
    self._constructFromIterable(iterable_simplices)

  def _constructFromIterable(self, iterable_simplices):
    for vertices in iterable_simplices:
      list(vertices).sort() # just in case
      simplex,_ = self.addSimplex(vertices)
      self._addFaces(simplex)

  def _addFaces(self,simplex): 
    if simplex.dim > 0:
      for removed_vertex in simplex.vertices:
        face,new_face = self.addSimplex([v for v in simplex.vertices if v != removed_vertex])
        simplex.faces.append(face) 
        face.cofaces.append(simplex)
        if new_face:
          self._addFaces(face)


  @classmethod
  def fromFile(cls, filename):
    ''' Constructs the complex from file

    Each row in the file corresponds to a simplex.
    It should be a comma-separated indexes of the 
    simplex. For example to enter [0,1,2] and [0,3],
    the file should only have 2 lines:

    0,1,2
    0,3

    The order in which the indecies appear does not
    matter. The order of the dimensions do not matter
    either (i.e. it is okay if 0,3 came first). 
   
    By definition, each subset of the simplex must also
    be in the complex, so we create all the implied 
    simplices. For example, the above file will create:
    [[0],[1],[2],[3],[0,1],[1,2],[0,2],[0,3],[1,2,3]]
    '''
    simplex_list = []
    with open(filename) as f:
      raw = f.read()
    for line in raw.split('\n'):
      if line != "":
        vertices = line.split(',')
        simplex_list.append(vertices)
    return cls(simplex_list)
  

  @classmethod
  def fromList(cls, simplex_list):
    ''' Constructs the complex from list of lists

    Each item in the list should contain an integer
    list with the indexes of the vertices in the simplex.
    Example:
       iterable_simplices = [[0,1,2],[0,3]]
    
    By definition, each subset of the simplex must also
    be in the complex, so we create all the implied 
    simplices. For example, [[0,1,2],[0,3]] will create:
    [[0],[1],[2],[3],[0,1],[1,2],[0,2],[0,3],[1,2,3]]
    '''
    return cls(simplex_list)


  def getNumberOfSimplexes(self):
    ''' Returns the total number of simplices.  '''
    s = 0
    for d_simplices in self.simplices:
      s+=len(d_simplices)
    return s

  def __str__(self):
    return str([d_simplex.values() for d_simplex in self.simplices])

  def __repr__(self):
    return str(self)
  

