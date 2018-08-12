# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

export vfields

# A type to make iteration over the simulation data more easy.
# 
# A Fields type can be constructed either with an array of times at
# which the simulation data is needed or without, in which case 
# iteration is performed over all possible snapshots, in `times(sim)`.
# 
# Type Parameters
# ---------------
# F : the type of the objects returned by this iterator. It can be a 
#     ScalarField or VectorField, and needs to be a type parameter so 
#     that we can define an eltype function and have efficient code. 
#     If this was not specified in the type parameter we would not 
#     know the return type and we would have bad performance.
struct FieldsIterator{F<:AbstractField,
                      S<:SimulationData,
                      V<:AbstractVector} <: AbstractVector{F}
      sim::S      # the simulation data object
    field::Symbol # a symbol for the field we want to load
       ts::V       # return snapshots at these times
    function FieldsIterator{F}(sim::S,
                            field::Symbol,
                                ts::V) where {S<:SimulationData,
                                              V<:AbstractVector{<:Real}, F}
        issubset(ts, times(sim)) || throw(SnapshotError("snapshots for all" *
                                    " times in `ts` must be available"))
        new{F, S, V}(sim, field, ts)
    end
end

# Constructor for vector fields iterator - default loads all velocity snapshots
vfields(sim::SimulationData{D, T},
      field::Symbol,
         ts::AbstractVector) where {D, T} = 
    FieldsIterator{VectorField{D, T, typeof(mesh(sim))}}(sim, field, ts)

# TODO
# Constructor for scalar fields iterator - default loads all scalar snapshots
# sfields(sim::SimulationData{D, T},
#          ts::AbstractVector=times(sim),
#         field::Symbol=:p) where {D}
#     FieldsIterator{ScalarField{D, T, typeof(mesh(sim))}}(field, sim, times(sim))

# Get the times at which snapshots will be given.
times(fs::FieldsIterator) = fs.ts

# Read only array interface
Base.size(fs::FieldsIterator) = (length(fs.ts), )
Base.getindex(fs::FieldsIterator{F}, i::Int) where {F} =
    _load_snapshot_vector(fs.sim, fs.ts[i], fs.field)::F

# TODO
# Base.getindex(fs::FieldsIterator{<:ScalarField}, i::Int) = 
    # _load_snapshot_scalar(fs.sim, fs.ts[i], fs.field)