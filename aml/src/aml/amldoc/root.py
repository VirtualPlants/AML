# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Root
	Root of the branching system containing a vertex
  Usage
	Root(v)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
	RestrictedTo (STRING): cf. Father
	ContainedIn (VTX): cf. Father
  Returned Object
	VTX
  Description
	This function is equivalent to Ancestors(v, EdgeType-> '<')@-1. It thus returns the root of the branching system containing the argument. This function never returns Undef.
  Background
	MTGs
  See Also
	
	MTG, 
	Extremities.

Examples
	AML> Ancestors(v) # set of ancestors of v
	<ARRAY(VTX)>[102,78,35,33,24,12]
	AML>  Root(v) # root of the branching system containing v
	<VTX>12


White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
"""
