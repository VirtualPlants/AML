# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
DateSample
	Array of observation dates of a vertex.
  Usage
	DateSample(v)
	DateSample(v, MinDate-> d1, MaxDate-> d2)
  Arguments
	v (VTX) : vertex of the active MTG.
  Optional Arguments
    MinDate (DATE) : defines a minimum date of interest.
    MaxDate (DATE) : defines a maximum date of interest.
  Returned Object
	ARRAY(DATE)
  Description
	Returns the set of dates at which a given vertex (passed as an argument) has been observed as an array of ordered dates. Options can be specified to define a temporal window and the total list of observation dates will be truncated according to the corresponding temporal window.
  Background
	Dynamic MTGs.
  See Also
	
	MTG, 
	FirstDefinedFeature, 
	LastDefinedFeature, 
	PreviousDate, 
	NextDate.
"""
