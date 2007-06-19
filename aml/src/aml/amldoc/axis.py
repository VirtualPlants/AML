# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Axis
	Array of vertices constituting a botanical axis
  Usage
	Axis(v)
      Axis(v, Scale-> s)
  Arguments
	v (VTX) : Vertex of the active MTG
  Optional Arguments
	Scale (STRING): scale at which the axis components are required.
  Returned Object
	ARRAY(VTX)
  Description
	An axis is a maximal sequence of vertices connected by `<'-type edges. Axis return the array of vertices representing the botanical axis which the argument v belongs to. The optional argument enables the user to choose the scale at which the axis decomposition is required.
  Background
	MTGs
  See Also
	
	MTG, 
	Path, 
	Ancestors.
  



White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
"""
