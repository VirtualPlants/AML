# This file has been extracted from AMAPmod documentation (refermanual 18)

# Predecessor
"""
    Father of a vertex connected to it by a '<' edge

    This function is equivalent to Father(v, EdgeType= '<'). 
    It thus returns the father (at the same scale) of the argument 
    if it is located in the same botanical. 
    If it does not exist, None is returned.

  
    Usage
    -----
    -----
        Predecessor(v)
    Parameters
    ----------
        - v (VTX) : vertex of the active MTG

    Optional Parameters
    -------------------
        - RestrictedTo (STRING): cf. Father
        - ContainedIn (VTX): cf. Father

    Returns
    -------
        vtx_id

  
    See Also
    --------
    --------
    
        `MTG`, `Father`, `Successor`.

    Examples
    --------
        >>> Predecessor(v)
        7
        >>> Father(v, EdgeType='+')
        >>> Father(v, EdgeType= '<')
        7
"""
