# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
Path
	Array of vertices defining the path between two vertices
  Usage
	Path(v1, v2)
  Arguments
	v1 (VTX) : vertex of the active MTG
	v2 (VTX) : vertex of the active MTG
  Returned Object
	ARRAY(VTX)
  Description
	This function returns the array of vertices defining the path between two vertices that are in an ancestor relationship. The vertex v1 must be an ancestor of vertex v2. Otherwise, if both vertices are valid, then the empty array is returned and if at least one vertex is undefined, Undef is returned.
  Note
	v1 can be equal to v2.
  Background
	MTGs
  See Also
	
	MTG, 
	Axis, 
	Ancestors.
  
Examples
	AML> v # print the value of v
	<VTX>78
	AML>  Ancestors(v)
	<ARRAY(VTX)>[78,45,32,10,4]
	AML>  Path(10,v)
	<ARRAY(VTX)>[10,32,45,78]
	AML>  Path(9,v) # 9 is not an ancestor of 78
	<ARRAY(ANY)>[]



White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
"""
