""" Stratification workload and the map of simplex -> strata

Distributed under MIT License

Authors: Ryo Asai (ryo@colfax-intl.com)
         Jay Shah (jshah3@nd.edu)
"""
import time
import sys
from CanonicalStratification.SimplicialComplex import *

class StrataMap:
  ''' Class containing map | simplex -> strata
  
  This class is the output of the stratify()
  function. Using getStratum() you can get the stratum 
  information of a SimplexNode object (see 
  SimplicialComplex.py). The function also accepts 
  a vertex set. For example:
    stratum = strata_map.getStratum([0,1,2]):

  Inidividal stratum currently only knows its index 
  and dimension. Index is simply given by an internal
  counter, and the dimension is the top dimension
  of the iteration when the stratum was created.

  The map is stored as a dict in simplex_stratum_map. 
  The dict key is the hash_val of simplex (see 
  SimplicialComplex.py)
  '''
  class Stratum:
    ''' Class representing Stratum.

    Do not manually create. Instead use createStratum()
    '''
    def __init__(self, index, top_dim):
      self.index = index
      self.top_dim   = top_dim


  def createStratum(self, top_dim):
    ''' Creates a new Stratum

    Increments an internal counter for the index.
    '''
    stratum = self.Stratum(self.strata_count, top_dim)
    self.strata_count+=1
    return stratum

  def getStratum(self, simplex):
    ''' Returns the stratum object for the given simplex. 

    If the stratum is not set for the simplex, it will return 
    None instead.

    Accepts SimplicialComplex.SimplexNode object or
    a list of vertices. For example:
      stratum = strata_map.getStratum([0,1,2]):
    '''
    if isinstance(simplex, SimplicialComplex.SimplexNode):
      simplex_hash = simplex.hash_val
    else:
      simplex_hash = SimplicialComplex.getSimplexHash(simplex)
    if simplex_hash in self.simplex_stratum_map:
      ret = self.simplex_stratum_map[simplex_hash]
    else:
      ret = None
    return ret


  def _setStratumOfSimplex(self, simplex, stratum):
    self.simplex_stratum_map[simplex.hash_val] = stratum

  def __init__(self):
    self.strata_count = 0
    self.simplex_stratum_map={}



def stratify(sc):
  ''' Main stratification algorithm

  The input is a SimplicialComplex object
  from SimplicialComplex.py. We require that
  the dimension is less than or equal to 3

  The algorithm itself is an iterative 
  procedure that is described in the associated
  publication (see README). Where applicable,
  the section that discusses the code will
  be referenced as "Ref {section number}".
  '''

  # Check on the dimension of the input
  assert len(sc.simplices)<=4, "Only complexes of dimension 3 or less is supported"

  def _getCofaces(sm, simplex, top_dim):
    ''' Returns immediate cofaces of a simplex
    that is in the current complex.

    Because we do not actually remove the simplices
    from the complex, we compute the current immediate
    cofaces on the fly.

    Ref 5.1
    '''
    cofaces = []
    for c in simplex.cofaces:
      s = sm.getStratum(c)
      # A coface is "current" it is either unassigned, or is assigned
      # to a stratum of the current top dimension.
      if (s is None) or (s.top_dim == top_dim):
        cofaces.append(c)
    return cofaces

  def _codimZeroOneCase(sm, c0_simplices, c1_simplices, top_dim):
    '''  The Codimension 0/1 case for checking stratum membership

    The workload boils down to a connected component analysis with
    an additional condition. 

    After the connected component analysis, we add all remaining 
    c_0 simplices to the current strata.

    Ref 5.3.1
    '''

    # Start connected component analysis for all unassigned simplices
    for simplex in c1_simplices:
      # We use stratum assignment for 'visited' check
      if sm.getStratum(simplex) is None:
        cofaces = _getCofaces(sm, simplex, top_dim)
        # Check if the simplex is the generic strata
        if len(cofaces) == 2:
          # Create a new strata for the connected component search
          stratum = sm.createStratum(top_dim)
          sm._setStratumOfSimplex(simplex, stratum)
          sm._setStratumOfSimplex(cofaces[0], stratum)
          sm._setStratumOfSimplex(cofaces[1], stratum)
          # Begin connected component search
          try:
            for neighbor in cofaces[0].faces:
              _connectedComponentSearch(sm, stratum, neighbor)
            for neighbor in cofaces[1].faces:
              _connectedComponentSearch(sm, stratum, neighbor)
          except RuntimeError: 
            # Python has a 999 recursion depth limit. It is possible
            #  to convert this to a while loop, but we are keeping it to recursion for now.
            raise RuntimeError('''Maximum recursion reached in connected component analysis. 
There is likely a large connected component with more than 999 elements. 
Increasing the default recursion limit of 999 with sys.setrecursionlimit(limit) is necessary for this complex.''')

    # Add all unassigned c0_simplex to the generic strata
    for simplex in c0_simplices:
      # We use stratum assignment for 'visited' check
      if sm.getStratum(simplex) is None:
        # New, isolated stratum for the simplex
        stratum = sm.createStratum(top_dim)
        sm._setStratumOfSimplex(simplex, stratum)

  def _connectedComponentSearch(sm, stratum, simplex):
    ''' Depth-first connected component analysis for the Codimension 0/1 
    case for checking stratum membership
    '''

     # we use stratum assignment for 'visited' check
    if sm.getStratum(simplex) is None:
      cofaces = _getCofaces(sm, simplex, stratum.top_dim)
      # Check if the simplex is the generic strata
      if len(cofaces) == 2:
        # add the simplex to the connected component
        sm._setStratumOfSimplex(simplex, stratum)
        sm._setStratumOfSimplex(cofaces[0], stratum)
        sm._setStratumOfSimplex(cofaces[1], stratum)
        # Begin connected component search
        for neighbor in cofaces[0].faces:
          _connectedComponentSearch(sm, stratum, neighbor)
        for neighbor in cofaces[1].faces:
          _connectedComponentSearch(sm, stratum, neighbor)
   
  def _codimTwoCase(sm, c2_simplices, top_dim):
    '''  The Codimension 2 case for checking stratum membership

    Each simplex is tested to see if all of its immediate cofaces
    are in the same stratum.

    Ref 5.3.2
    '''
    for simplex in c2_simplices:
      # We use stratum assignment for 'visited' check
      if sm.getStratum(simplex) is None:
        # New, isolated stratum for the simplex
        stratum = _uniqueStratumAmongCofaces(sm, simplex, top_dim)
        if stratum is not None:
          sm._setStratumOfSimplex(simplex, stratum)

  def _uniqueStratumAmongCofaces(sm, simplex, top_dim):
    ''' If all cofaces lie in the same stratum, return the stratum
    otherwise return None.

    Ref 5.2.1
    '''
    cofaces =  _getCofaces(sm, simplex, top_dim)
    stratum_set = set([sm.getStratum(c) for c in cofaces])
    if len(stratum_set) == 1:
      return stratum_set.pop()
    else:
      return None

  def _codimThreeCase(sm, c3_simplices, top_dim):
    '''  The Codimension 2 case for checking stratum membership

    Each simplex is first tested to see if all of its immediate 
    cofaces are in the same stratum. If they are, then the
    Euler Characteristic of the small link is computed.

    Ref 5.3.3
    '''
    for simplex in c3_simplices:
      # We use stratum assignment for 'visited' check
      if sm.getStratum(simplex) is None:
        # New, isolated stratum for the simplex
        stratum = _uniqueStratumAmongCofaces(sm, simplex, top_dim)
        if stratum is not None:
          sl = _getSmallLink(sm, simplex, top_dim)
          if (len(sl[0])-len(sl[1])+len(sl[2])) == 2:
            sm._setStratumOfSimplex(simplex, stratum)

  # Find the small link of the simplex
  def _getSmallLink(sm, simplex, top_dim):
    ''' Finding the small link of a simplex

    Workload is similar to that of connected component analysis
    in a directed graph.

    In the publication, there was an additional cutoff_dim
    parameter, but for codim <=3 the cutoff_dim is the same
    as codim. So it is omitted here.

    Ref 5.2.2
    '''

    # Small link array. Organized by dimension
    sl = [[] for i in range(top_dim-simplex.top_dim)]
    _addCofaceToSL(sl, sm, simplex, 0, top_dim)
    return sl

  def _addCofaceToSL(sl, sm, simplex, sl_dim, top_dim):
    for c in _getCofaces(sm, simplex, top_dim):
      if c not in sl[sl_dim]:
        sl[sl_dim].append(c)
        _addCofaceToSL(sl, sm, c, sl_dim+1, top_dim)


  ''' Main iterative part of the stratification. 

  Ref 5
  '''
  t0 = time.time()
  sm = StrataMap()

  for n_cur in reversed(range(len(sc.simplices))):
    # The forloop in the publication over k has been unrolled.

    #  Codim 0/1 case is called every time, but when dim is 0, c1_simplices is empty
    if n_cur-1 >= 0:
      _codimZeroOneCase(sm, sc.simplices[n_cur].values(), sc.simplices[n_cur-1].values(), n_cur)
    else:
      _codimZeroOneCase(sm, sc.simplices[n_cur].values(), [], n_cur)

    if n_cur-2 >= 0:
      _codimTwoCase(sm, sc.simplices[n_cur-2].values(), n_cur)

    if n_cur-3 >= 0:
      _codimThreeCase(sm, sc.simplices[n_cur-3].values(), n_cur)

  t1 = time.time()
  print("Stratification Time: %fs" % (t1-t0,))
  return sm


