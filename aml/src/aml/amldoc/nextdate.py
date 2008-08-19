# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
NextDate
        Next date at which a vertex has been observed after a specified date

    Usage
    -----
        NextDate(v, d)

    Parameters
    ----------
        v (VTX) : vertex of the active MTG.

    d (DATE) : departure date.

    Returns
    -------
        DATE

    Description
    -----------
        Returns the first observation date at which the vertex has been observed starting at date d and proceeding forward in time. None is returned if it does not exists.

    Background
    ----------
        Dynamic MTGs.

    See Also
    --------
        
        MTG, 
        DateSample, 
        FirstDefinedFeature, 
        LastDefinedFeature, 
        PreviousDate.
"""
