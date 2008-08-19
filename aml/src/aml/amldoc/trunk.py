# This file has been extracted from AMAPmod documentation (refermanual 18)

"""
Trunk
        Array of vertices constituting the bearing botanical axis of a branching system.

    Usage
    -----
        
                  Trunk(v)
      Trunk(v, Scale= s)
        

    Parameters
    ----------
        v (VTX) : Vertex of the active MTG.

    Optional Parameters
    -------------------
        Scale (STRING): scale at which the axis components are required.

    Returns
    -------
        ARRAY(VTX)

    Description
    -----------
        Trunk returns the array of vertices representing the botanical axis defined as the bearing axis of the whole branching system defined by v. The optional argument enables the user to choose the scale at which the trunk should be detailed.

    Background
    ----------
        MTGs

    See Also
    --------
        
        MTG, 
        Path, 
        Ancestors,
        Axis.
  



White vertices with double circle are argument(s) of the function

Black vertices are vertices returned by the function
 -=
"""
