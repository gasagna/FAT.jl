# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Simulations

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
       mesh

""" 
    OpenFoam simulation data object for post-processing.
    
    Parameters
    ----------
       casedir : the OpenFoam case directory;
     meshdtype : the type of the mesh points, face and cell 
                 centres, ... Defaults to Float64;
    fielddtype : the type of the field data. Defaults to Float64;

"""
type SimulationData{T<:Real}
    mesh::Mesh{T}
    fh::HDF5File
    t::Vector{Float64}
    dimensions::Int
    fielddtype::Type
end

function SimulationData(casedir::AbstractString; 
                        meshdtype::Type=Float64, 
                        fielddtype::Type=Float64)
    # Check that we are in a good openfoam directory
    iscasedir(casedir) || error("$casedir does not appear to " *
                                 "be an OpenFoam case directory")
    # check that the h5 file has been created in the casedir
    h5path = joinpath(casedir, "data.h5")
    ispath(h5path) || error("data.h5 not found in $(abspath(casedir)" *
                            "). Look at FAT.Utils.OpenFoamToHDF5 to " *
                            "create one.")
    # open the file
    fh = h5open(h5path)
    # read groups, which corresponds to snapshots
    t = sort!([parse(Float64, t) for t in names(fh)])
    # read mesh
    msh = Mesh(casedir, meshdtype)
    # Parse the number of dimensions of the simulation, either 2 or 3
    SimulationData(msh, 
                   fh, 
                   t, 
                   parse(Int, read(attrs(fh)["dimensions"])), 
                   fielddtype)
end

" Close hdf5 file "
close(sim::SimulationData) = close(sim.fh)

" Number of snapshots in the simulation "
length(sim::SimulationData) = length(sim.t)

" Get vector of times "
times(sim::SimulationData) = sim.t

" Number of spatial dimension of the fields "
ndims(sim::SimulationData) = sim.dimensions

" Get the mesh "
mesh(sim::SimulationData) = sim.mesh

" Type of data "
fielddtype(sim::SimulationData) = sim.fielddtype

" Nice printing "
function show(io::IO, sim::SimulationData)
    print(io, "OpenFOAM simulation object at $(object_id(sim)): \n")
    print(io, "  ~ data file $(abspath(filename(sim.fh)))\n")
    print(io, "  ~ $(length(sim.t)) snapshots available, from: ")
    print(io, times(sim)[1], " to ", times(sim)[end])
    print(io, "\n")
    print(io, "  ~ Mesh information:\n")
    print(io, "    ~ ")
    show(io, mesh(sim); space="     ")
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
    internalField = convert(Matrix{fielddtype(sim)}, 
                            read(sim.fh["$ts/W/internalField"]))
    boundaryField = convert(Matrix{fielddtype(sim)}, 
                            read(sim.fh["$ts/W/boundaryField"]))
    # Note that we take slices to avoid copying data. We could 
    # slice the h5 datasets directly, but there is bug #267 of 
    # HDF.jl that make slicing by column return a 
    # Matrix and not a Vector.
    VectorField(ntuple(i->ScalarField{ndims(sim), 
                                      fielddtype(sim), 
                                      typeof(mesh(sim))}(
                                     slice(internalField, :, i), 
                                     slice(boundaryField, :, i),
                                     mesh(sim)), ndims(sim)), mesh(sim))
end

function load_snapshot(sim::SimulationData, t::Real, var::Type{Val{:U}})
    ts = @sprintf "%.6f" t
    internalField = convert(Matrix{fielddtype(sim)}, 
                            read(sim.fh["$ts/U/internalField"]))
    boundaryField = convert(Matrix{fielddtype(sim)}, 
                            read(sim.fh["$ts/U/boundaryField"]))
    VectorField(ntuple(i->ScalarField{ndims(sim), 
                                      fielddtype(sim), 
                                      typeof(mesh(sim))}(
                                     slice(internalField, :, i), 
                                     slice(boundaryField, :, i),
                                     mesh(sim)), ndims(sim)), mesh(sim))
end

function load_snapshot(sim::SimulationData, t::Real, var::Type{Val{:ω}})
    ts = @sprintf "%.6f" t
    internalField = read(sim.fh["$ts/vorticity/internalField"])
    boundaryField = read(sim.fh["$ts/vorticity/boundaryField"])
    ScalarField{ndims(sim), 
                fielddtype(sim), 
                typeof(mesh(sim))}(internalField, 
                                   boundaryField,
                                   mesh(sim))
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