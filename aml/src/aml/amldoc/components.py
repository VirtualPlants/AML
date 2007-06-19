# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Components
	Set of components of a vertex
  Usage
	Components(v)
	Components(v, Scale-> 2)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
	Scale (INT) : scale of the components.
  Returned Object
	ARRAY(VTX)
  Description
	The set of components of a vertex is returned as an array of vertices. If s defines the scale of v, components are defined at scale s + 1. The array is empty if the vertex has no components. The order of the components in the array is not significant.
  Details
	When a scale is specified using optional argument Scale, it must be necessarily greater than the scale of the argument.
  Background
	MTGs
  See Also
	MTG, 
	Complex.



White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
"""
