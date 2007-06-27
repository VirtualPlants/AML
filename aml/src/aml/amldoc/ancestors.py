# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
Ancestors
	Array of all vertices which are ancestors of a given vertex
  Usage
	Ancestors(v)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
	RestrictedTo (STRING): cf. Father
      ContainedIn (VTX): cf. Father
      EdgeType (CHAR): cf. Father
  Returned Object
	ARRAY(VTX).
  Description
	This function returns the array of vertices which are located before the vertex passed as an argument. These vertices are defined at the same scale as v. The array starts by v, then contains the vertices on the path from v back to the root (in this order) and finishes by the tree root.
  Note
	The anscestor array always contains at least the argument vertex v.
  Background
	MTGs
  See Also
	
	MTG, 
	Descendants.

Examples
	AML> v # prints vertex v
	<VTX>78
	AML> Ancestors(v) # set of ancestors of v at the same scale
	<ARRAY(VTX)>[78,45,32,10,4]
	AML> Invert(Ancestors(v)) # To get the vertices in the 
	order from the root to the vertex v
	<ARRAY(VTX)>[4,10,32,45,78]
"""
