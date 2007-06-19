# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Coord
	Top coordinates of the geometric model of a component
  Usage
	Coord(p, v)
  Arguments
	p (PLANTFRAME) : plantframe containing the geometric representation of v.
	v (VTX) : vertex of the active MTG
  Returned Object
	ARRAY(REAL)
  Description
	This function returns the 3D-coordinates of the top of the box containing the geometric model of a plant component. The result is an array of 3 reals. 
  Note
	This function is similar to TopCoord. 
  Background
	MTGs
  See Also
	
	MTG, 
	PlantFrame, 
	TopCoord, 
	BottomCoord.
"""
