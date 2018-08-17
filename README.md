# CanonicalStratification

Library for canonical stratification of simplicial complexes.
This code is a companion to the publication "Algorithmic Stable Stratifications of Simplicial Complexes".

# Installation

There are no special module dependency for the core library.
However, "benchmark.py" requires `scipy` and `numpy` for producing Delaunay triangulation of 2-sphere and 3-ball.

To install:
	$ `pip install /path/to/CanonicalStratification`

# Usage

## TL;DR

There are two example files, `example_fromfile.py` and `example_fromlist.py` to try the stratification on a pre-made simplicial complex.
Also, `benchmark.py` will let you try the stratification on 2-sphere and 3-ball with an arbitrary number of vertices.

## More detailed explanation

You must first create a simplicial complex using the SimplicialComplex object.
This can be done through a python list or a file, and is discussed later.
But first we describe how a simplex is defined.

We assume that all vertices in the input dataset are indexed.
A simplex in the complex is defined by a set of indices for vertices.
For example, a simplex [0,1,2] consists of vertices indexed 0, 1 and 2. 

SimplicialComplex object is constructed by adding these simplices one by one.
For convenience, adding a simples adds all the _implied simplices_ as well.
Implied simplices of a simplex are all subsets of the simplex; mathematically they are the simplices in a closed sieve with the simplex at the top.
For example, a simplex [0,1,2] implies simplices \[0,1\],\[0,2\],\[1,2\],\[0\],\[1\] and \[2\]. 
These are added automatically. 

## Loading from lists

SimplicialComplex object can be loaded from an iterable of lists.
Each list in the iterable corresponds to a simplex, and contain the vertices of the simplex.
For example, the following specifies two simplicies \[0,1,2\] and \[0,3\]

`simplices = [[0,1,2],[0,3]]`

And a complex can be construced by.

`sc = SimplicialComplex.fromList(simplices)`

Note that this complex will have 9 simplices instead of 2 because of the implied simplices (see above).

Aditionally, the function takes any iterable of lists as an input. 
So, for example, you can input a generator of lists to this function as well.

## Loading from file

SimplicialComplex object can be loaded from a csv file.
Each line in the file corresponds to a simplex, and contain the vertices of the simplex.
For example, a file with the following entries specifies two simplicies \[0,1,2\] and \[0,3\]

`0,1,2`
`0,3`

And a complex can be construced by.

`sc = SimplicialComplex.fromFile(filename)`

Where `filename` is the string containing path to the file.

Note that this complex will have 9 simplices instead of 2 because of the implied simplices (see above).
