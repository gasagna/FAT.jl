# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

module FAT

# TODO
# ~ move the mesh data structure to the h5 file for faster loading

include("constants.jl")
include("ofio.jl")
include("meshes.jl")
# include("utils.jl")
include("field.jl")
# include("pod.jl")
include("simulations.jl")
# include("sparse.jl")
# include("galerkin.jl")

using .Constants, .Meshes, .Fields, .Simulation# .OFIO, .POD, .Sparse, .Galerkin


end
