# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
PreviousDate
        Previous date at which a vertex has been observed before a specified date

    Usage
    -----
        PreviousDate(v, d)

    Parameters
    ----------
        v (VTX) : vertex of the active MTG.
        d (DATE) : departure date.

    Returns
    -------
        DATE

    Description
    -----------
        Returns the first observation date at which the vertex has been observed starting at date d and proceeding backward in time. None is returned if it does not exists.

    Background
    ----------
        Dynamic MTGs.

    See Also
    --------
        
        MTG, 
        DateSample, 
        FirstDefinedFeature, 
        LastDefinedFeature, 
        NextDate.
"""
