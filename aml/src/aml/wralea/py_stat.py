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

#//////////////////////////////////////////////////////////////////////////////
# Input/output functions
#//////////////////////////////////////////////////////////////////////////////

def py_compound_ascii( filename ):
    if filename and filename != '':
        return (Compound(filename),)

#//////////////////////////////////////////////////////////////////////////////

def py_compound( sum_dist, dist ):
    if sum_dist and dist:
        return (Compound(sum_dist, dist),)

#//////////////////////////////////////////////////////////////////////////////

def py_convolution_ascii( filename ):
    if filename and filename != '':
        return (Convolution(filename),)

#//////////////////////////////////////////////////////////////////////////////

def py_convolution( dist ):
    if dist:
        return (Convolution(dist),)

#//////////////////////////////////////////////////////////////////////////////

def py_distribution_ascii( filename ):
    if filename and filename != '':
        return (Distribution(filename),)

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

def py_histogram( seq = [] ):
    """\
    Histogram([int]) -> Histogram
    Input:
        list of int
    Output:
        Histogram model
    """
    
    seq= self.get_input_by_key("sequence")
    if seq:
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



