# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Simulations

import DataStructures: SortedDict

import Base: getindex,
             start,
             next,
             done,
             length,
             show,
             ndims,
             eltype,
             size,
             linearindexing,
             LinearFast,
             close

import HDF5: h5read,
             HDF5File,
             h5open,
             attrs,
             filename

import FAT.Meshes: Mesh,
                   ncells,
                   nfaces,
                   nboundaryfaces

import FAT.Fields: AbstractField,
                   ScalarField,
                   VectorField,
                   mesh,
                   add!,
                   mul!

import FAT.OFIO: iscasedir

export SimulationData,
       fields,
       times,
       mesh,
       casedir

""" 
    OpenFoam simulation data object for post-processing.
    
    Parameters
    ----------
       casedir : the OpenFoam case directory
          mesh : the mesh object
          tmap : a dictionary of time/folder for the available snapshots

"""
type SimulationData{T<:Real}
    casedir::ASCIIString
    mesh::Mesh{T}
    tmap::SortedDict{Float64, ASCIIString, Base.Order.ForwardOrdering}
end

function SimulationData(casedir::AbstractString)
    # Check that we are in a good openfoam directory
    iscasedir(casedir) || error("$casedir does not appear to " *
                                 "be an OpenFoam case directory")
    
    # read snapshots and mesh
    tmap = get_folders_map(casedir)
    mesh = Mesh(casedir)

    SimulationData(casedir, mesh, tmap)
end

""" Get a dict of time=>folder with the times and 
    names of available snapshots. The keys of this
    dictionary, the snapshot times, are sorted in 
    ascending order.
"""
function get_folders_map(casedir::AbstractString)
    tmap = SortedDict{Float64, ASCIIString, Base.Order.ForwardOrdering}()
    for folder in readdir(casedir)
        # check that folder is a valid snapshot folder (parseable as a float)
        t = tryparse(Float64, folder)
        if !isnull(t) 
            tmap[get(t)] = folder
        end
    end
    tmap
end

" Get vector of times "
times(sim::SimulationData) = collect(keys(sim.tmap))

" Get the mesh "
mesh(sim::SimulationData) = sim.mesh

" Get the casedir "
casedir(sim::SimulationData) = sim.casedir

" Nice printing "
function show(io::IO, sim::SimulationData)
    print(io, "OpenFOAM simulation data object \n")
    print(io, "  ~ case directory $(abspath(casedir(sim)))\n")
    print(io, "  ~ $(length(sim.tmap)) snapshots available, from: ")
    print(io, times(sim)[1], " to ", times(sim)[end])
    print(io, "\n")
    print(io, "  ~ Mesh information:\n")
    print(io, "    ~ ")
    show(io, mesh(sim); gap="     ")
end

#= Load a velocity snapshot from file.

    Notes
    -----
    ~ We assume that it is a valid time, because this function is
    not mean to be used directly, and we fail badly. Use the 
    FieldsIteratorfor a general use.

    ~ We have a Val{:U} argument for type stability. If we have 
    other fields that we want to load, we define other load_snapshot 
    methods, with different Val{Symbol} arguments, so that every 
    method has a unique return type.
=#
function load_snapshot(sim::SimulationData, t::Real, var::Type{Val{:W}})
    ts = @sprintf "%.6f" t
    internalField = convert(Matrix{ftype(sim)}, 
                            read(sim.fh["$ts/W/internalField"]))
    boundaryField = convert(Matrix{ftype(sim)}, 
                            read(sim.fh["$ts/W/boundaryField"]))
    # Note that we take slices to avoid copying data. We could 
    # slice the h5 datasets directly, but there is bug #267 of 
    # HDF.jl that make slicing by column return a 
    # Matrix and not a Vector.
    VectorField(ntuple(i->ScalarField{ndims(sim), 
                                      ftype(sim), 
                                      typeof(mesh(sim))}(
                                      slice(internalField, :, i), 
                                      slice(boundaryField, :, i),
                                      mesh(sim)), ndims(sim)), mesh(sim))
end


# custom exception when snapshot is absent
type SnapshotError <: Exception 
    msg::AbstractString
end

function getindex(sim::SimulationData, t::Float64, var::Symbol) 
    t in sim.t || throw(SnapshotError("time $t not available"))
    return load_snapshot(sim, t, Val{var})
end

function getindex(sim::SimulationData, i::Int, var::Symbol)
    i > 0 && i <= length(sim.t) || throw(SnapshotError("invalid index $i"))
    return load_snapshot(sim, sim.t[i], Val{var})
end

# include FieldIterator
include("fieldsiterator.jl")

end