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

    nf = Factory( name= "Load_from_file", 
                  description= "Construction of a statistical object from an ASCII file.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "PyObjectFromFile",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Convolution", 
                  description= "Construction of an object of type convolution from an elementary distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_convolution",
                  )
    package.add_factory( nf )


    nf = Factory( name= "Binomial", 
                  description= "Construction of a binomial distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_binomial",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Poisson", 
                  description= "Construction of a poisson distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_poisson",
                  )
    package.add_factory( nf )

    nf = Factory( name= "NegativeBinomial", 
                  description= "Construction of a negative binomial distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_negativebinomial",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Uniform", 
                  description= "Construction of a uniform distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_dist_uniform",
                  )
    package.add_factory( nf )

    nf = Factory( name= "HiddenMarkov", 
                  description= "Construction of an object of type hidden_markov from an ASCII file.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_hiddenmarkov",
                  inputs=(dict(name='filename', interface=IFileStr),
                          dict(name='length', interface=IInt, value=20) ),
                  )
    package.add_factory( nf )

    nf = Factory( name= "HiddeniSemiMarkov", 
                  description= "Construction of an object of type hidden_semi-markov from an ASCII file.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_hiddensemimarkov",
                  inputs=(dict(name='filename', interface=IFileStr),
                          dict(name='length', interface=IInt, value=20),
                          dict(name='counting', interface=IBool, value=True)),
                  )
    package.add_factory( nf )

    nf = Factory( name= "Markov", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_markov",
                  inputs=(dict(name='filename', interface=IFileStr),
                          dict(name='length', interface=IInt, value=20) ),
                  )
    package.add_factory( nf )

    # TODO Mixture with a set of dist and proba

    '''
    '''

    nf = Factory( name= "Renewal", 
                  description= "Construction of a (either ordinary or equilibrium) renewal process from an inter-event distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "PyRenewal",
                  )
    package.add_factory( nf )

    nf = Factory( name= "Renewal_from_file", 
                  description= "Construction of a (either ordinary or equilibrium) renewal process from an inter-event distribution.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "PyRenewalAscii",
                  )
    package.add_factory( nf )

    nf = Factory( name= "SemiMarkov", 
                  description= "Construction of a semi-Markov chain from an ASCII file.", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_semimarkov",
                  inputs=(dict(name='filename', interface=IFileStr),
                          dict(name='length', interface=IInt, value=20),
                          dict(name='counting', interface=IBool, value=True)),
                  )
    package.add_factory( nf )
    nf = Factory( name= "Sequences", 
                  description= "", 
                  category = "STAT", 
                  nodemodule = "py_stat",
                  nodeclass = "py_sequences",
                  inputs=(dict(name='seq', interface=ISequence, value=[]),
                          dict(name="identifiers", interface=ISequence),
                          dict(name="indexParameter", 
                               interface=IEnumStr(["Position","Time"]), 
                               value = "Position")),
                  outputs=(dict(name="sequences", interface=ISequence),)
                  )
    package.add_factory( nf )

    '''
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
                  inputs=(dict(name='sequence',interface=ISequence,value=[]),),
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

