###############################################################################
# -*- python -*-
#
#       amlPy function implementation
#
#       Copyright or (C) or Copr. 2006 INRIA - CIRAD - INRA  
#
#       File author(s): Christophe Pradal <christophe.prada@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################

__doc__="""
amlPy functions
"""

__license__= "Cecill-C"
__revision__=" $Id$ "

#//////////////////////////////////////////////////////////////////////////////

from openalea.core import *
import openalea.aml as amlPy
#import amlPy
from amlPy import Class,Index,Scale,Order,Rank,Height


def max_mtg_scale(g):
    if not g:
        return 0
    
    g_str= str(g)
    scale= (g_str.split(',')[2]).strip()
    max_scale= int(scale.strip("levelnb="))
    return max_scale

#//////////////////////////////////////////////////////////////////////////////

class MTG(Node):
    """
MTG( filename, ErrorNb= 10, VtxNumber= 10000 ) -> MTG

Input: filename
Ouput : MTG object if the parsing process succeeds.
    """

    def __init__(self):

        Node.__init__(self)

        self.add_input( name = "filename", interface = IFileStr)
            
        self.add_output( name = "MTG", interface = None) 

        

    def __call__(self, inputs):
        """ inputs is the list of input values """

        from amlPy import MTG as mtg

        # We prefer here to get the value by key
        fn = self.get_input_by_key("filename")
        if fn:
            return (mtg(fn),)
        else:
            return (None,)

#//////////////////////////////////////////////////////////////////////////////

class VtxList(Node):
    """\
VtxList( mtg, Scale= 2 ) -> [vtx]

Input:
  MTG (MTG)
Optional parameters:
  Scale: used to select components at a particular scale.
Output:
  VtxList: vertex array
    """

    def __init__(self):

        Node.__init__(self)

        self.add_input( name = "MTG", interface = None) 
        self.add_input( name = "Scale", interface = IInt, value= 0)
            
        self.add_output( name = "VtxList", interface = None) 

        

    def __call__(self, inputs):
        """ inputs is the list of input values """

        from amlPy import VtxList as vertices
        from amlPy import Active, Activate

        # We prefer here to get the value by key
        g= self.get_input_by_key("MTG")
        if not g:
            g= Active()
        else:
            Activate(g)

        if not g:
            return ([],)
        
        scale = self.get_input_by_key("Scale")

        max_scale= max_mtg_scale(g)

        vtxs= []
        if scale == 0:
            vtxs= vertices()
        elif scale < max_scale:
            vtxs= vertices(Scale=scale)
        
        return (vtxs,)

#//////////////////////////////////////////////////////////////////////////////

class VtxFunction( Node ):
    """\
VtxFunction(string) -> func on vertex
Input:
  Name of the function.
Output:
  function that can be applied on a vertex.
    """
    
    vtx_func= { "Class" : Class,
                "Index" : Index,
                "Scale" : Scale,
                "Order" : Order,
                "Rank" : Rank,
                "Height" : Height }
    
    
    def __init__(self):
    
        Node.__init__(self)

        funs= self.vtx_func.keys()
        funs.sort()
        self.add_input( name = "Name", interface = IEnumStr(funs), value = funs[0]) 
        self.add_output( name = "VtxFunction", interface = None)

    def __call__(self, inputs):
        func_name= self.get_input_by_key("Name")
        return (self.vtx_func.get(func_name,None),)

#//////////////////////////////////////////////////////////////////////////////

class TopoFunction( Node ):
    """\
TopoFunction(string) -> func on vertex
Input:
  Name of the function.
Output:
  function that can be applied on a vertex.
    """
    
    vtx_func= { "Father" : amlPy.Father,
                "Successor" : amlPy.Successor,
                "Predecessor" : amlPy.Predecessor,
                "Root" : amlPy.Root,
                "Complex" : amlPy.Complex,
                "Location" : amlPy.Location,
                "Sons" : amlPy.Sons,
                "Ancestors" : amlPy.Ancestors, 
                "Descendants" : amlPy.Descendants,
                "Extremities" : amlPy.Extremities,
                "Components" : amlPy.Components,
                "ComponentRoots" : amlPy.ComponentRoots,
                "Axis" : amlPy.Axis,
                "Trunk" : amlPy.Trunk}
    
    
    def __init__(self):
    
        Node.__init__(self)

        funs= self.vtx_func.keys()
        funs.sort()
        self.add_input( name = "name", interface = IEnumStr(funs), value = funs[0]) 
        self.add_output( name = "f", interface = None)

    def __call__(self, inputs):
        func_name= self.get_input_by_key("name")
        return (self.vtx_func.get(func_name,None),)

#//////////////////////////////////////////////////////////////////////////////

class UnaryVtxFunc( Node ):
    
    def __init__(self,f):
    
        Node.__init__(self)
        self.f= f
        
        self.add_output( name = "f", interface = IFunction)

    def __call__(self, inputs):
        g= lambda x: self.f(x,*inputs)
        return (g,)

class Complex(UnaryVtxFunc):

    def __init__(self):
        UnaryVtxFunc.__init__(self,amlPy.Complex)
        self.add_input(name= "Scale", interface= IInt, value= 1 )

    def __call__(self, inputs):
        scale= self.get_input_by_key("Scale")
        g= lambda x: self.f(x,Scale= scale)
        return (g,)



class Components(UnaryVtxFunc):

    def __init__(self):
        UnaryVtxFunc.__init__(self,amlPy.Components)
        self.add_input(name= "Scale", interface= IInt, value= 1 )

    def __call__(self, inputs):
        scale= self.get_input_by_key("Scale")
        g= lambda x: self.f(x,Scale= scale)
        return (g,)


#//////////////////////////////////////////////////////////////////////////////

class Plot( Node ):
    """\
Plot(aml object) -> Plot the object
Input:
  aml object
  
    """
    
    def __init__(self):
    
        Node.__init__(self)

        self.add_input( name = "obj", interface = None, value = []) 

    def __call__(self, inputs):
        obj= self.get_input_by_key("obj")
        if obj:
          amlPy.Plot(obj)



