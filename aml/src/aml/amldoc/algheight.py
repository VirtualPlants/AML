# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
AlgHeight
	Algebraic value defining the number of components between two components
  Usage
	AlgHeight(v1, v2)
  Arguments
	v1 (VTX) : vertex of the active MTG
      v2 (VTX) : vertex of the active MTG
  Returned Object
	INT. If v1 is not an ancestor of v2 (or vice versa) , or if v1 and v2 are not defined at the same scale, an error value Undef is returned.
  Description
	This function is similar to function Height(v1, v2): it returns the number of components between two components, at the same scale, but takes into account the order of vertices v1 and v2. The result is positive if v1 is an ancestor of v2, and negative if v2 is an ancestor of v1.
  Background
	MTGs
  See Also
	
	MTG, 
	Rank, 
	Order, 
	Height, 
	EdgeType, 
	AlgOrder,
	AlgRank.
"""
