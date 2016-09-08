# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
Simulate
	Generation of a random sample from a distribution.
  Usage
	Simulate(dist, size)
	Simulate(mixt, size)
	Simulate(convol, size)
	Simulate(compound, size)
	
  Arguments
	dist (distribution),
	mixt (mixture),
	convol (convolution),
	compound (compound),
	size (int): sample size.
	
  Returned Object
	
	 If the first argument is of type distribution and if 0 < size � 1000000, an object of type HISTOGRAM is returned, otherwise no object is returned.
	 If the first argument is of type mixture and if 0 < size � 1000000, an object of type mixture_data is returned, otherwise no object is returned.
	 If the first argument is of type convolution and if 0 < size � 1000000, an object of type convolution_data is returned, otherwise no object is returned.
	 If the first argument is of type compound and if 0 < size � 1000000, an object of type compound_data is returned, otherwise no object is returned.
	The returned object of type HISTOGRAM, mixture_data, convolution_data or compound_data contains both the simulated sample and the model used for simulation.
	
  See Also
	
	Distribution,
	Mixture,
	Convolution,
	Compound,
	ExtractHistogram.
"""
