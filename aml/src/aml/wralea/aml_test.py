
mtg10 = MTG(amapmod_dir + "databases/MTGFiles/AppleTree/wij10.mtg")

gus = VtxList(Scale->3)

#############################################
#
#  Exploration functions
#
#############################################

# 1. sampling function
# number of internodes per entity

innb(_x) = Size(Components(_x, Scale->4))

# 2. Sample
# set of growth units for a given order 

gus_o(_order) = Foreach _x In gus : Select(_x, Order(_x) == _order)

# Iteration : construction of the sample of values

sample_innb_gu(_order) = Foreach _x In gus_o(_order) : innb(_x)

histo_nb_gu0 = Histogram(sample_innb_gu(0))
Plot(histo_nb_gu0)
