###############################################################################
# -*- python -*-
#
#       amlPy function definition 
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
vplants.stat function definition
"""

__license__= "Cecill-C"
__revision__=" $Id$ "

from openalea.core.external import *

def define_factory(package):
    """ Define factories for arithmetics nodes """

    # Input/output functions
    nf = Factory( name= "Compound", 
                  description= "Construction of an object of type compound from a sum distribution and an elementary distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_compound",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Compound_from_file", 
                  description= "Construction of an object of type compound from an ASCII file.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_compound_ascii",
                  inputs=tuple(dict(name='filename',interface=IFileStr)),
                  )
    package.add_factory( nf )

    nf = Factory( name= "Convolution", 
                  description= "Construction of an object of type convolution from an elementary distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_convolution",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Convolution_from_file", 
                  description= "Construction of an object of type convolution from an ASCII file.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_convolution_ascii",
                  inputs=tuple(dict(name='filename',interface=IFileStr)),
                  )
    package.add_factory( nf )

    nf = Factory( name= "Distribution", 
                  description= "Construction of a parametric discrete distribution (either binomial, Poisson, negative binomial or uniform) from a file.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_distribution_ascii",
                  inputs=tuple(dict(name='filename',interface=IFileStr)),
                  )
    package.add_factory( nf )

    nf = Factory( name= "Distribution_Binomial", 
                  description= "Construction of a binomial distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_binomial",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Distribution_Poisson", 
                  description= "Construction of a poisson distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_poisson",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Distribution_NegativeBinomial", 
                  description= "Construction of a negative binomial distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_negativebinomial",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Distribution_Uniform", 
                  description= "Construction of a uniform distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_uniform",
                  )
    package.add_factory( nf )

    '''
    nf = Factory( name= "HiddenMarkov", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_hidenmarkov",
                  )
    package.add_factory( nf )

    nf = Factory( name= "HiddeniSemiMarkov", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_hidensemimarkov",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Markov", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_markov",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Mixture", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_mixture",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Renewal", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_renewal",
                  )
    package.add_factory( nf )

    nf = Factory( name= "SemiMarkov", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_semimarkov",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Sequences", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_sequences",
                  )
    package.add_factory( nf )

    nf = Factory( name= "TimeEvents", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_timeevents",
                  )
    package.add_factory( nf )

    nf = Factory( name= "TopParameters", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_topparameters",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Tops", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_tops",
                  )
    package.add_factory( nf )

    nf = Factory( name= "VectorDistance", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_vectordistance",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Vectors", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_vectors",
                  )
    package.add_factory( nf )

    nf = Factory( name= "", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_",
                  )
    package.add_factory( nf )

    '''
    nf = Factory( name= "Histogram", 
                  description= "Histogram from a sequence of int.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_histogram",
                  inputs=tuple(dict(name='sequence',interface=ISequence,value=[])),
                  outputs=tuple(dict(name='histogram')),
                  )
    package.add_factory( nf )

    nf = Factory( name= "Plot", 
                  description= "Graphical output of an object of the STAT module using the GNUPLOT software.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_plot",
                  outputs=(),
                  )
    package.add_factory( nf )

