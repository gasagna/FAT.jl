# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

""" A type to make iteration over the simulation data more easy. Look 
    at the test_Simulations.jl for examples on how to use this thing.

    A Fields type can be constructed either with an array of times at
    which the simulation data is needed or without, in which case 
    iteration is performed over all possible snapshots, in `times(sim)`.

    Type Parameters
    ---------------
    F : the type of the objects returned by this iterator. It can be a 
        ScalarField or VectorField, and needs to be a type parameter so 
        that we can define an eltype function and have efficient code. 
        If this was not specified in the type parameter we would not 
        know the return type and we would have bad performance.
"""
type FieldsIterator{F<:AbstractField} <: AbstractVector{F}
    sim::SimulationData
    var::Symbol
    ts::Vector{Float64}
    function FieldsIterator(sim::SimulationData,
                            var::Symbol,
                            ts::AbstractVector)
        issubset(ts, sim.t) || throw(SnapshotError("snapshots for all" *
                                    " times in `ts` must be available"))
        new(sim, var, Float64[ts;])
    end
end

" Get the times at which snapshots will be given. "
times(fs::FieldsIterator) = fs.ts

# implement array interface for FieldsIterator type, see 
# http://docs.julialang.org/en/release-0.4/manual/
# interfaces/#abstract-arrays
getindex(fs::FieldsIterator, i::Int) = 
    load_snapshot(fs.sim, fs.ts[i], Val{fs.var})
linearindexing(::Type{FieldsIterator}) = LinearFast()
size(fs::FieldsIterator) = length(fs.ts,)
length(fs::FieldsIterator) = length(fs.ts)
eltype{F}(fs::FieldsIterator{F}) = F

for (fname, fftype) in zip([:(:U), :(:W)], 
        [:VectorField, :VectorField])
    @eval begin 
        fields(sim::SimulationData, ::Type{Val{$fname}}) = 
            FieldsIterator{$fftype{ndims(sim), 
                                 ftype(sim), 
                                 typeof(mesh(sim))}}(sim, $fname, sim.t)
        fields(sim::SimulationData, 
                  ::Type{Val{$fname}}, 
                ts::AbstractVector) = 
            FieldsIterator{$fftype{ndims(sim), 
                                  ftype(sim), 
                                  typeof(mesh(sim))}}(sim, $fname, ts)
    end
end