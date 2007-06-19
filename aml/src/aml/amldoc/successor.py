# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Successor
	Vertex that is connected to a given vertex by a '<' edge type (i.e. in the same botanical axis).
  Usage
	Successor(v)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
	RestrictedTo (STRING): cf. Father
	ContainedIn (VTX): cf. Father
  Returned Object
	VTX
  Description
	This function is equivalent to Sons(v, EdgeType-> '<')@1. It returns the vertex that is connected to a given vertex by a `<' edge type (i.e. in the same botanical axis). If many such vertices exist, an arbitrary one is returned by the function. If no such vertex exists, Undef is returned.
  Background
	MTGs
  See Also
	
	MTG, 
	Sons, 
	Predecessor.

Examples
	AML> Sons(v) # 
	<ARRAY(VTX)>[3,45,47,78,102]
	AML>  Sons(v, EdgeType-> '+') # set of vertices borne by v
	<ARRAY(VTX)>[3,45,47,102]
	AML>  Sons(v, EdgeType-> '<') # set of successors of v
	<ARRAY(VTX)>[78]
	AML>  Successor(v)
	<VTX>78
"""
