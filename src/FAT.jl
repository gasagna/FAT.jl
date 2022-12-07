# ------------------------------------------------------------------- #
# Copyright 2015-2019, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

module FAT

# You need to have these modules
# - DataStructures
# - HeterogeneousVectors

include("constants.jl")
include("ofio.jl")
include("meshes.jl")
include("field.jl")
include("fieldsiterator.jl")
include("operators.jl")
# include("pod.jl")
include("simulations.jl")
# include("galerkin.jl")

end
