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
# Adapters
def adapt2list(arg):
    if arg is None:
        return []
    elif type(arg) == types.ListType:
        return arg
    else:
        return [arg]

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
    l=adapt2list(list_of_dist)
    if l:
        return Convolution(*l)

def py_mixture( list_of_dist ):
    l=adapt2list(list_of_dist)
    if l:
        arg=[]
        for w,d in l:
            arg.append(w)
            arg.append(d)

    return Mixture(*arg)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_binomial( inf_bound= 0, sup_bound = 10, proba = 0.5 ):
    return (Distribution("BINOMIAL",inf_bound, sup_bound, proba),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_poisson( inf_bound= 0, param = 10 ):
    return (Distribution("P",inf_bound, param),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_negativebinomial( inf_bound= 0, param = 10., proba = 0.5 ):
    return (Distribution("NB",inf_bound, param, proba),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_uniform( inf_bound= 0, sup_bound = 10 ):
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
        if type(obj) == types.ListType:
            Plot(*obj)
        else:
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

#//////////////////////////////////////////////////////////////////////////////
# Compare family : TODO

def py_compare_frequency( histos, Type ):
    _types = {"Numeric" : "N","Ordinal" : "O", "Symbolic" : "S"}
    _type=_types[Type]
    if histos is None:
        return None,
    args = list(adapt2list(histos))
    args.append(_type)
    
    return Compare(*args),
        

def py_compare_vectors( vectors, vector_distance ):
    if vectors is None or vector_distance is None:
        return None,
    return Compare(vectors, vector_distance),

#//////////////////////////////////////////////////////////////////////////////

def py_comparisontest(Type, histo1, histo2):
    if Type in ['F', 'T', 'W'] and histo1 and histo2:
        return ComparisonTest(Type, histo1, histo2)

#//////////////////////////////////////////////////////////////////////////////

def py_estimate_dist( histo,
                      distribution,
                      MinInfBound,
                      InfBoundStatus,
                      ):
    if not histo: return
    if distribution == 'NON-PARAMETRIC':
        return Estimate(histo, distribution)
    else:
        return Estimate(histo, 
                        distribution, 
                        MinInfBound=MinInfBound, 
                        InfBoundStatus=InfBoundStatus)


def py_estimate_mixture( histo,
                         components,
                         MinInfBound,
                         InfBoundStatus,
                         DistInfBoundStatus,
                         NbComponents,
                         Penalty,):

    if not histo: return
    args = [histo]
    args.append("MIXTURE")

    mix = components.split(',')
    mix = map(lambda x:x.strip(),mix)
    args.extend(mix)

    kwds = {}
    kwds['MinInfBound'] = MinInfBound
    kwds['InfBoundStatus'] = InfBoundStatus
    kwds['DistInfBoundStatus'] = DistInfBoundStatus
    kwds['NbComponent'] = NbComponents

    if NbComponents == "Estimated":
        kwds['Penalty'] = Penalty

    return Estimate(*args, **kwds)

def py_estimate_conv( histo,
                      dist,
                      Estimator,
                      InitialDistribution,
                      MinInfBound,
                      NbIteration,
                      Penalty,
                      Weight,
                      Outside,
                    ):

    if not histo or not dist: return
    
    kwds = {}
    kwds['Estimator']=Estimator
    kwds['MinInfBound']=MinInfBound
    kwds['NbIteration']=NbIteration
    if InitialDistribution:
        kwds['InitialDistribution'] = InitialDistribution

    if Estimator == "PenalizedLikelihood":
        kwds['Penalty'] = Penalty
        kwds['Weight'] = Weight
        kwds['Outside'] = Outside
    
    return Estimate(histo,dist,*kwds)
"""
def py_estimate_distrib( histo,
                      distribution,
                      mixtures,
                      unknow,
                      MinInfBound,
                      InfBoundStatus,
                      DistInfBoundStatus,
                      NbComponents,
                      Penalty,
                      Parametric,
                      InitialDistribution ):
    if not histo: return
    args = [histo]
    args.append(distribution)

    if distribution == 'MIXTURE':
        mix = mixtures.split(',')
        mix = map(lambda x:x.strip(),mix)
        args.extend(mix)

    if distribution in ['CONVOLUTION','COMPOUND'] and unknow != '': 
        args.append(unknow)

    kwds = {}
    if MinInfBound != 0:
        kwds['MinInfBound'] = MinInfBound

    if distribution != 'NON-PARAMETRIC':
        if InfBoundStatus != 'Free':
            kwds['InfBoundStatus'] = InfBoundStatus

    if distribution == 'MIXTURE':
        if DistInfBoundStatus != 'Fixed':
            kwds['DistInfBoundStatus'] = DistInfBoundStatus
        if NbComponents == "Estimated":
            kwds['NbComponent'] = NbComponents
            kwds['Penalty'] = Penalty

    if distribution in ['CONVOLUTION','COMPOUND']:
        kwds['Parametric'] = Parametric
        if InitialDistribution:
            kwds['InitialDistribution'] = InitialDistribution
            if 'MinInfBound' in kwds:
                del kwds['MinInfBound']

    return Estimate(*args, **kwds)
"""

def py_merge( data ):
    
    if data:
        if type(data) == types.ListType:
            return Merge(*data)
        else:
            return Merge(data)

def py_extractdata( model ):
    if model is not None:
        return ExtractData(model)

def py_shift(obj, param=0):
    if obj: return Shift(obj,param)

def py_shiftn(obj, variable=0, param=0):
    if obj: return Shift(obj,variable, param)

def py_cluster(obj, mode, step, information_ratio, limits, AddVariable):
    if obj is None: return
    if mode == "Step":
        return Cluster(obj, mode, step)
    elif mode == 'Information':
        return Cluster(obj, mode, information_ratio)
    elif mode == "Limit":
        _limits = adapt2list(limits)
        return Cluster(obj, mode, _limits)

def py_simulate_dist(obj, size=100):
    if obj is None:
        return
    return Simulate(obj,size)
