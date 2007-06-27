# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
ComponentRoots
	Set of roots of the tree graphs that compose a vertex
  Usage
	ComponentRoots(v)
	ComponentRoots(v, Scale-> s)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
	Scale (STRING): scale of the component roots.
  Returned Object
	ARRAY(VTX)
  Description
	In a MTG, a vertex may have be decomposed into components. Some of these components are connected to each other, while other are not. In the most general case, the components of a vertex are organized into several tree-graphs. This is for example the case of a MTG containing the description of several plants: the MTG root vertex can be decomposed into tree graphs (not connected) that represent the different plants. This function returns the set of roots of these tree graphs at scale Scale(v)+1. The order of these roots is not significant.
  Details
	When a scale different from Scale(v)+1 is specified using the optional argument Scale, this scale must be greater than that of the vertex argument.
  Background
	MTGs
  See Also
	
	MTG, 
	Components,
	Trunk.

Examples
	>>> v=MTGRoot() # global MTG root
	<ARRAY(VTX)>0
	>>> ComponentRoots(v) # set of first vertices at scale 1
	<ARRAY(VTX)>[1,34,76,100,199,255]
	>>> ComponentRoots(v, Scale> 2) # set of first vertices at scale 2
	<ARRAY(VTX)>[2,35,77,101,200,256]



White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
"""
