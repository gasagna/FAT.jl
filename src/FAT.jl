isdefined(Base, :__precompile__) && __precompile__(true)

module FAT

# TODO
# ~ make the Field and SimulationData aware of the format the simulation is stored 
#   in the hdf file, by making them parametric
# ~ move the mesh data structure to the h5 file for faster loading
# ~ when you open a simulation you should have a command to close it
# ~ add integrals over patches 

include("constants.jl")
include("ofio.jl")
include("meshes.jl")
include("utils.jl")
include("field.jl")
include("pod.jl")
include("simulations.jl")
include("sparse.jl")
include("galerkin.jl")

using .Constants, .Meshes, .Utils, .Fields, .Simulations, .OFIO, .POD, .Sparse
using .Galerkin


end
