# ------------------------------------------------------------------- #
# Copyright 2015-2022, Davide Lasagna, AFM, University of Southampton #
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
include("simulations.jl")

# unused, but still here
# include("pod.jl")
# include("galerkin.jl")

end
