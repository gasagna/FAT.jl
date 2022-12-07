 # ------------------------------------------------------------------- #
# Copyright 2015-2019, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
import DataStructures: SortedDict

export SimulationData,
       times,
       mesh,
       casedir,
       load_vector_snapshot

# ============================================================================ #
# SIMULATION DATA OBJECT
# ============================================================================ #
struct SimulationData{D,
                      T<:Real,
                      S<:SortedDict{Float64, String, Base.Order.ForwardOrdering}}
    casedir::String         # path of simulation data
       mesh::Mesh{T}        # mesh object
       dims::NTuple{D, Int} # load velocity along these directions
       tmap::S              # snapshot times
end

Base.ndims(s::SimulationData{D}) where {D} = D
Base.eltype(s::SimulationData{D, T}) where {D, T} = T

# set of monotonically increasing dimensions
const _ALLOWED_DIMENSIONS = [(1, 2, 3), (1, 2), (1, 3), (2, 3)]

# constructor
function SimulationData(casedir::AbstractString,
                           dims::NTuple{D, Int},
                               ::Type{T} = Float64) where {D, T}
    # check dimensions make sense
    dims in _ALLOWED_DIMENSIONS || 
        throw(ArgumentError("dimensions $dims not understood"))

    # Check that we are in a good openfoam directory
    iscasedir(casedir) || error("$casedir does not appear to " *
                                 "be an OpenFoam case directory")
    
    # read snapshots and mesh
    tmap = _get_folders_map(casedir)
    mesh = Mesh(casedir, T)

    SimulationData{D, T, typeof(tmap)}(casedir, mesh, dims, tmap)
end

" Get vector of times "
times(sim::SimulationData) = collect(keys(sim.tmap))

" Get the mesh "
mesh(sim::SimulationData) = sim.mesh

" Get the casedir "
casedir(sim::SimulationData) = sim.casedir

" Nice printing "
function Base.show(io::IO, sim::SimulationData)
    print(io, "OpenFOAM simulation data object \n")
    print(io, "  ~ case directory $(abspath(casedir(sim)))\n")
    print(io, "  ~ $(length(sim.tmap)) snapshots available, from: ")
    print(io, times(sim)[1], " to ", times(sim)[end])
    print(io, "\n")
    print(io, "  ~ Mesh information:\n")
    print(io, "    ~ ")
    show(io, mesh(sim); gap="     ")
    return nothing
end


# ============================================================================ #
# CUSTOM EXCEPTION WHEN SNAPSHOT IS NOT AVAILABLE
# ============================================================================ #
struct SnapshotError <: Exception 
    msg::AbstractString
end


# ============================================================================ #
# UTILITY FUNCTIONS FOR LOADING THE DATA OR OTHER STUFF
# ============================================================================ #

# Get a dictionary of time=>folder with times and names of available snapshots. 
# The keys (the snapshot times), are sorted in  ascending order.
function _get_folders_map(casedir::AbstractString)
    tmap = SortedDict{Float64, String}(Base.Order.Forward)
    for folder in readdir(casedir)
        # check that folder is a valid snapshot folder (parseable as a float)
        t = tryparse(Float64, folder)
        if t != nothing
            tmap[t] = folder
        end
    end
    return tmap
end


function load_vector_snapshot(sim::SimulationData{D}, 
                                t::Real, 
                              var::Symbol) where {D}
    # return a valid file, or raise an error
    filename = _get_filename(sim, t, var)
    # read boundary and internal fields
    intf, bndf = read_vector_field(sim.casedir, filename, sim.dims)
    # type parameters of the vector field
    params = length(sim.dims), Float64,  Mesh{Float64}
    # tuple of scalar fields, each containing one component
    scalars = ntuple(i->ScalarField(mesh(sim), D, intf[i], bndf[i]), D)
    # output
    return VectorField(mesh(sim), scalars)
end

# TODO
# function _load_snapshot_scalar(sim::SimulationData, t::Real, var::Symbol)
# end 

# Build filename to be read
function _get_filename(sim::SimulationData, t::Real, var::Symbol)
    t in keys(sim.tmap) || throw(SnapshotError("time $t not available"))
    filename = joinpath(sim.casedir, sim.tmap[t], string(var))
    ispath(filename) || throw(SnapshotError("field $var not available at time $t"))
    return filename
end