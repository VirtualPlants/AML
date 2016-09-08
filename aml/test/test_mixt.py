from openalea.aml import *

#########################################################################
#########################################################################
#
#  Frequency distributions
#
#  Objective: Analyzing the number of nodes of growth units in selected architectural
#             position considering the respective roles of preformation and neoformation,
#
#  Methods: comparison tests, one-way variance analysis,
#           estimation of finite mixture of distributions.
#
#  Wild cherry tree: number of nodes per growth unit (GU)
#
#  Data: Dominique Fournier
#
#  meri1.his: order 1,
#  meri1.his: order 2,
#  meri1.his: order 3, GU 1,
#  meri1.his: order 3, GU 2,
#  meri5.his: short shoots.
#
#
#  Poplar: number of nodes per growth unit
#
#  Data: Yves Caraglio and Herve Rey
#
#  peup1.his: order 2,
#  peup2.his: order 3,
#  peup3.his: order 4,
#  peup4.his: order 5,
#  peup5.his: order 3, GU 4,
#  peup6.his: order 3, acrotony.
#
#########################################################################
def test_mixt():
    meri1 = Histogram("data/meri1.his")
    meri2 = Histogram("data/meri2.his")
    meri3 = Histogram("data/meri3.his")
    meri4 = Histogram("data/meri4.his")
    meri5 = Histogram("data/meri5.his")
    
# Plot(meri1, meri2, meri3, meri4, meri5)

    Compare(meri1, meri2, meri3, meri4, meri5, "NUMERIC")
    
    ComparisonTest("F", meri1, meri2)
    ComparisonTest("T", meri1, meri2)
    ComparisonTest("W", meri1, meri2)
    
    ComparisonTest("F", meri1, meri3)
    ComparisonTest("T", meri1, meri3)
    ComparisonTest("W", meri1, meri3)
    
# estimation of a mixture of two distributions assuming a first sub-population of GUs
# made only of a preformed part and a second sub-population made of both a preformed part
# and a neoformed part

    mixt1 = Estimate(meri2, "MIXTURE", "B", "B")
    
    meri = Merge(meri1, meri2, meri3, meri4, meri5)
    
    # model selection approach: estimation of both the mixture parameters and
    # the number of components 

    mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",  NbComponent="Estimated")
    # mixt2 = Estimate(meri, "MIXTURE", "NB", "NB")
    # Plot(ExtractDistribution(mixt2, "Mixture"))
    # Display(mixt2)
    
    mixt_data = ExtractData(mixt2)
    
    
    dist5 = Estimate(meri5, "BINOMIAL")
    # Display(dist5, Detail=2)
    # Plot(dist5)
    
    histo5 = Simulate(dist5, 100)
    # Display(histo5, Detail=2)
    # Plot(histo5)
    
    
    peup1 = Histogram("data/peup1.his")
    peup2 = Histogram("data/peup2.his")
    peup3 = Histogram("data/peup3.his")
    peup4 = Histogram("data/peup4.his")
    peup5 = Histogram("data/peup5.his")
    peup6 = Histogram("data/peup6.his")
    
    mixt10 = Estimate(peup2, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")
    
    peup = Merge(peup1, peup2, peup3, peup4, peup5, peup6)
    
    histo1 = Shift(peup, -1)
    histo2 = Cluster(peup, "Information", 0.8)
    histo3 = Cluster(peup, "Step", 10)
    histo4 = Cluster(peup, "Limit", [13, 24])
    # Display(histo4, Detail=2)
    # Plot(histo4)
    

    mixt11 = Estimate(peup, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")
    # mixt11 = Estimate(peup, "MIXTURE", "B", "NB")
    
    assert mixt11
