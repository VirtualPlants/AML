# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
Descendants
	Set of vertices in the branching system borne by a vertex.
  Usage
	Descendants(v)
  Arguments
	v (VTX) : vertex of the active MTG
  Optional Arguments
    RestrictedTo (STRING): cf. Father
    ContainedIn (VTX): cf. Father
    EdgeType (CHAR): cf. Father
  Returned Object
	ARRAY(VTX)
  Description
	This function returns the set of descendants of its argument as an array of vertices. The array thus consists of all the vertices, at the same scale as v, that belong to the branching system starting at v. The order of the vertices in the array is not significant. 
  Note
	The argument always belongs to the set of its descendants.
  Background
	MTGs
  See Also
	
	MTG, 
	Ancestors.

Examples
    AML> v # 
    <VTX>78
    AML> Sons(v) # set of sons of v
    <ARRAY(VTX)>[78,99,101]
    AML> Descendants(v) # set of descendants of v
    <ARRAY(VTX)>[78,99,101,121,133,135,156,171,190]



White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
"""
