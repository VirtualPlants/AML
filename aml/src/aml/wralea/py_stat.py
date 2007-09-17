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
import types

#//////////////////////////////////////////////////////////////////////////////
# Input/output functions
#//////////////////////////////////////////////////////////////////////////////

def py_compound_ascii( filename ):
    if filename and filename != '':
        return (Compound(filename),)
def py_convolution_ascii( filename ):
    if filename and filename != '':
        return (Convolution(filename),)
def py_distribution_ascii( filename ):
    if filename and filename != '':
        return (Distribution(filename),)
def py_histogram_ascii(filename):
    if filename and filename != '':
        return (Histogram(filename),)


#//////////////////////////////////////////////////////////////////////////////

def py_compound( sum_dist, dist ):
    if sum_dist and dist:
        return (Compound(sum_dist, dist),)

#//////////////////////////////////////////////////////////////////////////////

def py_convolution( list_of_dist ):
    if hasattr(list_of_dist,'__iter__'):
        return (Convolution(list_of_dist),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_binomial( inf_bound= 2, sup_bound = 5, proba = 0.8 ):
    return (Distribution("BINOMIAL",inf_bound, sup_bound, proba),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_poisson( inf_bound= 0, param = 12.2 ):
    return (Distribution("P",inf_bound, param),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_negativebinomial( inf_bound= 5, param = 3.2, proba = 0.4 ):
    return (Distribution("NB",inf_bound, param, proba),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_uniform( inf_bound= 2, sup_bound = 5 ):
    return (Distribution("U",inf_bound, sup_bound),)

#//////////////////////////////////////////////////////////////////////////////

def py_hiddenmarkov( filename, length=20 ):

    if filename and filename != '':
        return (HiddenMarkov(filename,length),)

#//////////////////////////////////////////////////////////////////////////////

def py_hiddensemimarkov( filename, length=20, counting=True ):

    if filename and filename != '':
        return (HiddenSemiMarkov(filename,length,counting),)

#//////////////////////////////////////////////////////////////////////////////

def py_markov( filename, length=20 ):

    if filename and filename != '':
        return (Markov(filename,length),)

#//////////////////////////////////////////////////////////////////////////////

def py_mixture( filename ):

    if filename and filename != '':
        return (Mixture(filename),)

#//////////////////////////////////////////////////////////////////////////////

class PyRenewal(Node):
    def __init__(self):
        Node.__init__(self)
        self.Types=["Equilibriun","Ordinary"]
        self.add_input(name="inter_event")
        self.add_input(name="Type", interface = IEnumStr(self.Types), value = self.Types[0])
        self.add_input(name="ObservationTime", interface=IInt, value=20)

    def __call__(self, inputs):
        inter_event= self.get_input("inter_event")
        Type=self.get_input("Type")
        obs_time=self.get_input("ObservationTime")
        return (Renewal(inter_event,Type,obs_time),)

class PyRenewalAscii(Node):
    def __init__(self):
        Node.__init__(self)
        self.Types=["Equilibriun","Ordinary"]
        self.add_input(name="filename", interface=IFileStr)
        self.add_input(name="Type", interface = IEnumStr(self.Types), value = self.Types[0])
        self.add_input(name="ObservationTime", interface=IInt, value=20)

    def __call__(self, inputs):
        filename= self.get_input("filename")
        Type=self.get_input("Type")
        obs_time=self.get_input("ObservationTime")
        if filename:
            return (Renewal(filename,Type,obs_time),)

#//////////////////////////////////////////////////////////////////////////////

def py_semimarkov( filename, length=20, counting=True ):

    if filename and filename != '':
        return (SemiMarkov(filename,length,counting),)

#//////////////////////////////////////////////////////////////////////////////

def py_histogram( seq = [] ):
    """\
    Histogram([int]) -> Histogram
    Input:
        list of int
    Output:
        Histogram model
    """
    
    if hasattr(seq,'__iter__') and (len(seq)>0):
        return (Histogram(seq),)

#//////////////////////////////////////////////////////////////////////////////

def py_plot( obj ):
    """
    Plot(aml object) -> Plot the object
    Input:
        aml object
    """
    
    if obj:
        Plot(obj)

#//////////////////////////////////////////////////////////////////////////////

class PyObjectFromFile(Node):
    """ 
    Construct a model (Distribution, vectors, ...) from a file.
    """

    Types= {"BINOMIAL" : Distribution,
            "POISSON" : Distribution,
            "NEGATIVE_BINOMIAL" : Distribution,
            "UNIFORM" : Distribution,
            "COMPOUND" : Compound,
            "CONVOLUTION" : Convolution,
            "MIXTURE" : Mixture,
            "HISTOGRAM" : Histogram,
            "SEQUENCES" : Sequences,
            "TIME_EVENTS" : TimeEvents,
            "TOPS" : Tops,
            "VECTOR_DISTANCE" : VectorDistance,
            "VECTORS" : Vectors,
            
            }
    def __init__(self):
        Node.__init__(self)
        self.add_input(name="filename", interface=IFileStr)
        self.add_input(name="Type", 
                       interface = IEnumStr(self.Types.keys()), 
                       value = "BINOMIAL")
        self.add_output(name='Model')
        self.set_caption("import model")

    def __call__(self, inputs):
        filename= self.get_input("filename")
        name=self.get_input("Type")
        self.set_caption(name)

        klass = self.Types.get(name)
        if filename:
            return (klass(filename),)

#//////////////////////////////////////////////////////////////////////////////

def py_sequences(seq=[], identifiers=[], indexParameter='Position'):
    if seq:
        if identifiers:
            return (Sequences(seq,identifier,indexParameter),)
        else:
            return(Sequences(seq),)



