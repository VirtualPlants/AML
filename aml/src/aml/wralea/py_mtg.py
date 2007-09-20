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
from openalea.aml import *


def max_mtg_scale(g):
    if not g:
        return 0
    
    g_str= str(g)
    scale= (g_str.split(',')[2]).strip()
    max_scale= int(scale.strip("levelnb="))
    return max_scale

def set_mtg(g):
    old_g = Active()
    if not g:
        g= old_g
    elif g != old_g:
        Activate(g)

    return g

#//////////////////////////////////////////////////////////////////////////////

class py_MTG(Node):
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

        # We prefer here to get the value by key
        fn = self.get_input("filename")
        g = None
        if fn:
            try:
                g = MTG(fn)
            except Exception, e:
                print e
                g = None

        return (g,)

#//////////////////////////////////////////////////////////////////////////////

class py_VtxList(Node):
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

        #from openalea.aml import VtxList as vertices
        #from openalea.aml import Active, Activate

        # We prefer here to get the value by key
        g= self.get_input("MTG")
        g = set_mtg(g)
        if not g:
            return ([],)
        
        scale = self.get_input("Scale")

        max_scale= max_mtg_scale(g)

        vtxs= []
        if scale == 0:
            vtxs= VtxList()
        elif scale <= max_scale:
            vtxs= VtxList(Scale=scale)
        
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
        self.add_input( name = "Vtx" ) 
        self.add_output( name = "VtxFunction", interface = None)

    def __call__(self, inputs):
        func_name= self.get_input("Name")
        vtx = self.get_input("Vtx")
        f = self.vtx_func.get(func_name,None)
        self.set_caption(func_name)

        if not vtx:
            return f
        elif callable(vtx):
            return lambda x: f(vtx(x))
        else:
            return f(vtx)

#//////////////////////////////////////////////////////////////////////////////

class TopoFunction( Node ):
    """\
TopoFunction(string) -> func on vertex
Input:
  Name of the function.
Output:
  function that can be applied on a vertex.
    """
    
    vtx_func= { "Father" : Father,
                "Successor" : Successor,
                "Predecessor" : Predecessor,
                "Root" : Root,
                "Complex" : Complex,
                "Location" : Location,
                "Sons" : Sons,
                "Ancestors" : Ancestors, 
                "Descendants" : Descendants,
                "Extremities" : Extremities,
                "Components" : Components,
                "ComponentRoots" : ComponentRoots,
                "Axis" : Axis,
                "Trunk" :Trunk}
    
    
    def __init__(self):
    
        Node.__init__(self)

        funs= self.vtx_func.keys()
        funs.sort()
        self.add_input( name = "name", interface = IEnumStr(funs), value = funs[0]) 
        self.add_input( name = "Vtx" ) 
        self.add_output( name = "f", interface = None)

    def __call__(self, inputs):
        func_name= self.get_input("name")
        vtx = self.get_input("Vtx")
        f = self.vtx_func.get(func_name,None)
        self.set_caption(func_name)
        if not vtx:
            return f
        elif callable(vtx):
            return lambda x: f(vtx(x))
        else:
            return f(vtx)

#//////////////////////////////////////////////////////////////////////////////

class UnaryVtxFunc( Node ):
    
    def __init__(self,f):
    
        Node.__init__(self)
        self.f= f
        
        self.add_output( name = "f", interface = IFunction)

    def __call__(self, inputs):
        g= lambda x: self.f(x,*inputs)
        return (g,)

class py_Complex(UnaryVtxFunc):

    def __init__(self):
        UnaryVtxFunc.__init__(self,Complex)
        self.add_input(name= "Vtx", interface= IInt )
        self.add_input(name= "Scale", interface= IInt, value= 0 )

    def __call__(self, inputs):
        vtx = self.get_input("Vtx")
        scale= self.get_input("Scale")
        if vtx:
            return (self.f(vtx,Scale=scale),)
        else:
            g= lambda x: self.f(x,Scale= scale)
            return (g,)



class py_Components(UnaryVtxFunc):

    def __init__(self):
        UnaryVtxFunc.__init__(self,Components)
        self.add_input(name= "Vtx", interface= IInt )
        self.add_input(name= "Scale", interface= IInt, value= 0 )

    def __call__(self, inputs):
        vtx = self.get_input("Vtx")
        scale= self.get_input("Scale")
        if vtx:
            return (self.f(vtx,Scale=scale),)
        else:
            g= lambda x: self.f(x,Scale= scale)
            return (g,)


#//////////////////////////////////////////////////////////////////////////////

class py_Plot( Node ):
    """\
Plot(aml object) -> Plot the object
Input:
  aml object
  
    """
    
    def __init__(self):
    
        Node.__init__(self)

        self.add_input( name = "obj", interface = None, value = []) 

    def __call__(self, inputs):
        obj= self.get_input("obj")
        if obj:
          Plot(obj)

#//////////////////////////////////////////////////////////////////////////////

class py_PlantFrame( Node ):
    """\
Plot(aml object) -> Plot the object
Input:
  aml object
  
    """
    
    def __init__(self):
    
        Node.__init__(self)

        self.add_input( name = "MTG", interface = None) 
        self.add_input( name = "Vertex", interface = IInt, value= 0)
        self.add_input( name = "Scale", interface = IInt, value= 0)
        self.add_input( name = "DressingData")
        self.add_output(name = "plantframe")

    def __call__(self, inputs):
        g= self.get_input("MTG")
        vtx= self.get_input("Vertex")
        scale= self.get_input("Scale")
        d=self.get_input("DressingData")

        g = set_mtg(g)
        if not g:
            return (None,)
        
        max_scale= max_mtg_scale(g)
        if scale > max_scale:
            scale = max_scale
        
        kwds={}
        if scale:
            kwds["Scale"] = scale
        if d:
            kwds["DressingData"] = d

        try:
            pf = PlantFrame(vtx, **kwds)
        except Exception, e:
            print e
            pf=None

        return (pf,)

#//////////////////////////////////////////////////////////////////////////////

def py_dressingdata(g,filename):
    g =set_mtg(g)
    if not g:
        return

    if filename:
        try:
            d = DressingData(filename)
        except Exception, e:
            print e
            d = None
        return (d,)


