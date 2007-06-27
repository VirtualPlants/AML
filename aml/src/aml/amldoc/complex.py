# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
Complex
	Complex of a vertex
  Usage
	Complex(v)
	Complex(v, Scale-> 2)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
    Scale (INT) : scale of the complex.
  Returned Object
	VTX
  Description
	Returns the complex of the argument. The complex of a vertex v is has a scale lower than v: Scale(v)-1. In a MTG, every vertex except for the MTG root (cf. MTGRoot()), has a uniq complex. Undef is returned if the argument is the MTG root or if the vertex is undefined.
  Details
	When a scale different from Scale(v)-1 is specified using the optional argument Scale, this scale must be lower than that of the vertex argument.
  Background
	MTGs
  See Also
	
	MTG, 
	Components.
"""
