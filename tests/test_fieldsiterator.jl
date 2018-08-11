using Base.Test

using FAT.Fields
using FAT.Meshes
using FAT.Simulation

# do tests on ascii data format
sim = SimulationData("ldc_test_ascii", (1, 2))

# construct iterator object
us = vfields(sim)

# check lenngth
@test length(us) == 11