# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Predecessor
	Father of a vertex connected to it by a '<' edge
  Usage
	Predecessor(v)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
	RestrictedTo (STRING): cf. Father
	ContainedIn (VTX): cf. Father
  Returned Object
	VTX
  Description
	This function is equivalent to Father(v, EdgeType-> '<'). It thus returns the father (at the same scale) of the argument if it is located in the same botanical. If it does not exist, Undef is returned.
  Background
	MTGs
  See Also
	
	MTG, 
	Father, 
	Successor.

Examples
	AML> Predecessor(v)
	<VTX> 7
	AML>  Father(v, EdgeType-> '+')
	Undef
	AML>  Father(v, EdgeType-> '<')
	<VTX> 7
"""
