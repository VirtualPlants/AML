# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Sons
	Set of vertices borne or preceded by a vertex
  Usage
	Sons(v)
	Sons(v, EdgeType-> '+')
	Sons(v, Scale-> 3)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
	RestrictedTo (STRING): cf. Father
	ContainedIn (VTX): cf. Father
	EdgeType (CHAR) : filter on the type of sons.
	Scale (INT) : set the scale at which sons are considered.
  Returned Object
	ARRAY(VTX)
  Description
	The set of sons of a given vertex is returned as an array of vertices. The order of the vertices in the array is not significant. The array can be empty if there are no son vertices.
  Details
	When the option EdgeType is applied, the function returns the set of sons that are connected to the argument with the specified type of relation. Note that Sons(v, EdgeType-> '<') is not equivalent to Successor(v). The first function returns an array of vertices while the second function returns a vertex.
    	
	The returned vertices have the same scale as the argument. However, coarser or finer vertices can be obtained by specifying the the optional argument Scale at which the sons are considered.
  Background
	MTGs
  See Also
	
	MTG, 
	Father, 
	Successor, 
	Descendants.

Examples
	AML> Sons(v) # 
	<ARRAY(VTX)>[3,45,47,78,102]
	AML>  Sons(v, EdgeType-> '+') # set of vertices borne by v 
	<ARRAY(VTX)>[3,45,47,102]
	AML>  Sons(v, EdgeType-> '<') # set of successors of v on the same axis
	<ARRAY(VTX)>[78]



White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
"""
