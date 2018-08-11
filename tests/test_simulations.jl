using Base.Test

using FAT.Fields
using FAT.Meshes
using FAT.Simulation

casedirname = "ldc_test_ascii"

# use a test run
sim = SimulationData(casedirname, (1, 2))

# check times
@test times(sim) == collect(0:10)

# check we return a mesh
@test typeof(mesh(sim)) == Mesh{Float64}

# check casedirname
@test casedir(sim) == casedirname