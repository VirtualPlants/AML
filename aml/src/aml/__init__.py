# -#-python-#-
#  ----------------------------------------------------------------------------
#
#  amlPy package.
#
#  File author(s): Christophe Pradal, Cyril Bonnard
#
#  $Id$
#
#  ----------------------------------------------------------------------------
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  ----------------------------------------------------------------------------
#
#
#  $Id$
#
#  ----------------------------------------------------------------------------
#

# Version
from amlPy_version import amlPy_version

__version__= amlPy_version

#from __future__ import nested_scopes

try:
   import PlantGL 
except ImportError:
   pass

import _amlPymodule

_amlPy = _amlPymodule

# If we cannot import Numeric, it should have been detected at
# installation time and amlPy should have been compiled properly. So,
# don't complain.
try:
    import Numeric
except ImportError:
    pass


try:
    import openalea.plantgl.all
except ImportError:
   try:
      import PlantGL
   except ImportError:
      pass



###########################################################################
############################# Global names ################################
###########################################################################


###########################################################################
########################### Conversion system #############################
###########################################################################

   
# Manage the global mode

def setmode(mode):
    _amlPy.set_mode(mode)

def getmode():
    return _amlPy.get_mode()

setmode(1)

class AML(object):
    """ main class """

    # called when the attribute is not found,
    # eg always except for TRUE and FALSE
    def __getattr__(self, name):
        return self.__getitem__(name)

    # return the amobj from r.__dict__ at key=name
    def __getitem__(self, name):
        #print "__getitem__", name
        obj = self.__dict__[name] = self.__dict__.get(name, _amlPy.get(name))
        return obj

    # enables the instance to be callable eg: aml("AML code")
    def __call__(self, s):
        #print "__call__", s
        return self.eval(self.parse(text=s))

# main instance
aml = AML()

# hack: usual AMAPmod functions

# Equivalent AMAPmod function in Python
import math, os

Acos = math.acos
Asin = math.asin
Atan = math.atan
Cos = math.cos
Log = math.log
Log10 = math.log10
Sin = math.sin
Sqrt = math.sqrt
Tan = math.tan

Abs = abs

Pwd = os.getcwd
Cd = os.chdir
Dir = os.getcwd
#Ls=os.listdir(Dir())
Rm = os.remove
Mkdir = os.mkdir

Max = max
Min = min

def Swich(x, dict, default):
    """
    Switch(x, dict, default)
    Return dict[x] or default
    """
    return dict.get( x, default )


__liste = """
AMAPMOD_DIR
Abs
Activate
Active
AlgHeight
AlgOrder
AlgRank
AllPos 
Alpha  
Ancestors
Angle    
Append   
Array
Axis          
Azimuth       
Beta          
Black
Blue 
BottomCoord   
BottomDiameter
Ceil          
Class         
ClassScale    
Complex       
ComponentRoots
Components    
Coord         
Copy          
Cov
CurrentMTG
CurrentPlottedObj
CurrentWindow 
CWD
Date          
DateFormat
DateSample    
DateUnit
Day
DefaultWindow
Defined       
Delete        
Descendants   
Display       
DressingData  
EDist         
Echo          
EchoOff       
EchoOn        
EdgeType      
ElectricalModel
EulerAngles   
Exists        
Exp
ExpSequence   
Extract       
Extremities   
Father        
Feature       
FileHelp      
Filter        
FirstDefinedFeature
Flatten       
Floor         
GeomScene     
Green
Head          
Height        
HOME
Hour
Identity      
IfThenElse    
Index         
InsertAt      
Inter         
Interval      
Invert        
LastDefinedFeature
Length        
LightBlue
List          
Location
Ls
MSVoxel       
MTDModel      
MTG
MTGRoot       
MatchingExtract
Mean          
Minus         
Minute
Mod
NewPlot       
NextDate      
Norm          
Octree        
Order         
PDir          
Path          
Pi
PlantFrame    
Plot
Plus
Pos
Predecessor   
PreviousDate  
Print         
Prod          
ProdSeries    
Rank          
Red
RelBottomCoord
RelTopCoord   
Rand          
Remove        
RemoveAt      
Rename        
Rint          
Root
ROOTDIR
SDir          
SProd         
Save          
SaveMTG       
Scale         
Second
Select        
Series        
Set
SetAmapmodDir 
SetMinus      
SetViewer     
Size          
Sons          
Sort          
StartViewer   
StopViewer    
SubArray      
Successor     
Sum
Switch        
Tail          
Times         
ToArray       
ToInt         
ToList        
ToReal        
ToSet         
ToString      
ToTimeUnit    
TopCoord      
TopDiameter   
TreeMatching  
TreeMatchingLoad
Truncate      
Trunk         
Undef
Union         
VProd         
Var
ViewerHelp    
ViewerIsRunning
ViewerSelection
Violet
VirtualPattern
VtxList       
WaitViewer    
White
Yellow

AddAbsorbingRun
Cluster
Clustering
Compare
ComparisonTest
Compound
ComputeAutoCorrelation
ComputeCorrelation
ComputeInitialRun
ComputePartialAutoCorrelation
ComputeRankCorrelation
ComputeSelfTransition
ComputeStateSequences
ComputeWhiteNoiseCorrelation
ContingencyTable
Convolution
Cross
Cumulate
Difference
Distribution
Estimate
ExtractData
ExtractDistribution
ExtractHistogram
ExtractVectors
Fit
HiddenMarkov
HiddenSemiMarkov
HiddenVariableOrderMarkov
Histogram
IndexExtract
LengthSelect
Markov
Merge
MergeVariable
Mixture
ModelSelectionTest
MovingAverage
NbEventSelect
RecurrenceTimeSequences
Regression
RemoveApicalInternodes
RemoveRun
Renewal
RenewalData
Reverse
Segmentation
SegmentationExtract
SelectIndividual
SelectVariable
SemiMarkov
Sequences
Shift
Simulate
Symmetrize
Thresholding
TimeEvents
TimeScaling
TimeSelect
ToDistanceMatrix
ToDistribution
ToHistogram
TopParameters
Tops
Transcode
TransformPosition
TransitionCount
TruncateDistribution
Unnormalize
ValueSelect
VariableOrderMarkov
VariableScaling
VarianceAnalysis
VectorDistance
Vectors
WordCount
"""

__liste = __liste.split()
__dictionnary = globals()


for fct in __liste:
    __dictionnary[fct] = aml.__getattr__(fct)


def amlhelp(name):
   """ Display the AML documentation of the specified Function """

   try:
      baseimp = "openalea.aml.amldoc." + name.lower()
      docmodule = __import__(baseimp, globals(), locals(), ['__doc__'])
      help(docmodule)

   except Exception, e:
      print "Documentation unavailable:", e



for fct in __liste:
   baseimp = "openalea.aml.amldoc." + fct.lower()
   try:
      docmodule = __import__(baseimp, globals(), locals(), ['__doc__'])
   except:
      print fct, "NO DOC"
   
