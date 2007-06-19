# This file has been extracted from AMAPmod documentation (refermanual 18)

__doc__ = """
Sequences
	Construction of a set of sequences from multidimensional arrays of integers, from data generated by a renewal process or from an ASCII file.
  Usage
	Sequences(array1, Identifiers->[1, 8, 12])
	Sequences(arrayn, Identifiers->[1, 8, 12],
	IndexParameter->"Position")
	
	Sequences(timev)
	
	Sequences(file_name)	
  Arguments
	array1 (array(array(int))): input data for univariate sequences
	arrayn (array(array(array(int)))): input data for multivariate sequences,
	
	timev (renewal_data),
	
	file_name (string).
	
  Optional Arguments 
	Identifiers (array(int)): explicit identifiers of sequences. This optional argument can only be used if the first argument is of type array(array(int/array(int))).
	IndexParameter (string): type of the explicit index parameter: "Position" or "Time" (the default: implicit discrete index parameter starting at 0). This optional argument can only be used if the first argument is of type array(array(int/array(int))).
	
  Returned Object
	If the construction succeeds, an object of type sequences or discrete_sequences is returned, otherwise no object is returned. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is � 15.
  Description
	The data structure of type array(array(array(int))) should be constituted at the most internal level of arrays of constant size. If the optional argument IndexParameter is set at "Position" or "Time", the data structure of type array(array(array(int))) is constituted at the most internal level of arrays of size 1+n (index parameter, n variables attached to the explicit index parameter). If the optional argument IndexParameter is set at "Position", only the index parameter of the last array of size 1+n is considered and the first component of successive elementary arrays (representing the index parameter) should be increasing. If the optional argument IndexParameter is set at "Time", the first component of successive elementary arrays should be strictly increasing.
	
  See Also
	Save, 
	ExtractHistogram, 
	ExtractVectors, 
	AddAbsorbingRun,
	Cluster, 
	Cumulate, 
	Difference, 
	Indexextract, 
	Lengthselect, 
	Merge, 
	MergeVariable, 
	MovingAverage, 
	RecurrenceTimeSequences, 
	RemoveRun, 
	Reverse, 
	SegmentationExtract, 
	SelectIndividual, 
	SelectVariable, 
	Shift, 
	Transcode, 
	ValueSelect,
	VariableScaling, 
	ComputeCorrelation, 
	ComputePartialAutoCorrelation, 
	ComputeSelfTransition, 
	Compare (sequences), 
	Compare (Markovian models of seuqences), 
	Compare (Markovian models), 
	Estimate (Markovian models), 
	ComputeStateSequences, 
	Simulate (Markovian models).
"""
