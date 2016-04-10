# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Meshes

import DataStructures: DefaultDict

import Base: start,
             next,
             done,
             length,
             eltype,
             show

import FAT.OFIO: reader,
                 read_boundary,
                 iscasedir

export Patch, 
       patchname,
       firstfaceID,
       lastfaceID,
       nfaces

export Mesh,
       patch,
       npatches,
       patches,
       ncells,
       nfaces,
       nboundaryfaces,
       ninternalfaces,
       patchfaces,
       internalfaces,
       cellcentres

include("point.jl")
include("polygonalface.jl")
include("polyhedralcell.jl")
include("patches.jl")

""" Calculate interpolation weights from cell centres to face centre 

    Notes
    -----
    See formula 6.31 in the book The Finite Volume Method in 
    Computational Fluid Dynamics: An Advanced Introduction with 
    OpenFOAM® and Matlab. See figure 6.24

    Parameters
    ----------
    fcentre : face centre    
    ocentre : owner cell centre   
    ncentre : neighbour cell centre   
    fsvec   : face surface vector  
"""
function interpolationWeight{T<:Real}(fcentre::Point{T}, 
                                      ocentre::Point{T}, 
                                      ncentre::Point{T}, 
                                      fsvec::Point{T})
    down = (fcentre - ocentre)*fsvec
    dnei = (ncentre - fcentre)*fsvec 
    down/(down+dnei)
end

# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Mesh definition ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~
""" Type for meshes with arbitrary polyhedral cells 

    Type parameters
    ---------------
    T : type of the data, i.e. points, faces surface vectors, volumes, ...

    Fields
    ------
    fsvecs : face surface vectors
    fowners : integer indices representing the owner cell of each face
    fneighs : integer indices representing the neighbour cell of each face
    cvolumes : cell volumes
    ccentres : cell mass centroid
    patches : dictionary of patches
    αs : interpolation coefficients for each face
    nboundaryfaces : total number of boundary faces, including 
                     those on empty patches
    ninternalfaces : the total number of internal faces
"""
type Mesh{T}
    fsvecs::Vector{Point{T}}
    fowners::Vector{UInt32}
    fneighs::Vector{UInt32}
    cvolumes::Vector{T}
    ccentres::Vector{Point{T}}
    patches::Dict{Symbol, Patch}
    αs::Vector{T}
    nboundaryfaces::Int
    ninternalfaces::Int
    function Mesh(fcentres, fsvecs, fowners, fneighs, 
                  cvolumes, ccentres, patches)
        nboundaryfaces = sum([nfaces(v) for (p, v) in patches])
        ninternalfaces = length(fcentres) - nboundaryfaces
        αs = Vector{T}(ninternalfaces)
        for i = 1:ninternalfaces
            αs[i] = interpolationWeight(fcentres[i], 
                                        ccentres[fowners[i]], 
                                        ccentres[fneighs[i]], fsvecs[i])
        end
        new(fsvecs, fowners, fneighs, cvolumes, ccentres, patches, αs, 
            nboundaryfaces, ninternalfaces)
    end
end

# data type of the mesh points, areas, volumes, ...
eltype{T}(::Mesh{T}) = T

" Number of faces on the domain boundary "
@inline nboundaryfaces(m::Mesh) = m.nboundaryfaces

" Number of internal faces "
@inline ninternalfaces(m::Mesh) = m.ninternalfaces

" Returns a range object to iterate over all internal faces "
@inline internalfaces(m::Mesh) = 1:m.ninternalfaces

""" 
    Returns two range objects in a zip to iterate over the faces of a patch. 
    The first contains indices of the faces, so that we can index the 
    vectors `fowners` and `fneighs` appropriately. Note that these two 
    have length equal to the total number of faces in the mesh, i.e. 
    both the internal and external faces.  The second range object
    is used to index the vector `.boundaryField` in `ScalarField`
    object, which contains the solution at the centres of the boundary 
    faces. An example of usage of this function is the function 
    `der!` in src/field.jl.
"""
function patchfaces(m::Mesh, p::Patch) 
    faceIDs = firstfaceID(p):lastfaceID(p)
    return zip(faceIDs, faceIDs - ninternalfaces(m))
end

" Total number of faces in the mesh "
nfaces(m::Mesh) = length(m.fsvecs)

" Total number of cells in the mesh "
ncells(m::Mesh) = length(m.cvolumes)

" Dictionary of mesh patches "
patches(m::Mesh) = m.patches

" Get a patch named `pname`"
patch(m::Mesh, pname::Symbol) = m.patches[pname]

" Get `dir` coordinates of cell centres "
function cellcentres{T}(m::Mesh{T}, dir::Symbol)
    out = Vector{T}(ncells(m))
    cc = m.ccentres
    for i in 1:length(out)
        out[i] = getfield(cc[i], dir)
    end
    out
end

" Get cell centres "
cellcentres(m::Mesh) = 
    (cellcentres(m, :x), cellcentres(m, :y), cellcentres(m, :z))

function show(io::IO, mesh::Mesh; gap::AbstractString=" ")
    print(io, "Mesh object at $(object_id(mesh)):  \n")
    print(io, "$gap ~ $(ncells(mesh)) cells                 \n")
    print(io, "$gap ~ $(nfaces(mesh)) total faces           \n")
    print(io, "$gap ~ $(nboundaryfaces(mesh)) boundary faces\n")
    print(io, "$gap ~ $(length(keys(patches(mesh)))) patches: \n")
    for (patchname, v) in mesh.patches
        print(io, "$gap   ~ $(string(patchname)): $(nfaces(v)) faces\n")
    end
end

""" Read the mesh in an OpenFoam case directory 

    Parameters
    ----------
    casedir : the case directory
    dtype   : main datatype for everything; defaults to Float64

"""
function Mesh{T<:Real}(casedir::AbstractString, dtype::Type{T}=Float64)

    # check before loading
    iscasedir(casedir) || error("$casedir is not an OpenFoam case!")

    # read all data. These are 1-based, and not 0-based as the OF files
    points_data = reader(casedir, "points"; dtype=dtype)
    faces_data  = reader(casedir, "faces")
    fowners     = reader(casedir, "owner")
    fneighs     = reader(casedir, "neighbour")

    # size of mesh
    local nfaces = length(faces_data)
    local ncells = maximum(fowners)

    # we also need to construct
    fsvecs   = Vector{Point{dtype}}(nfaces)
    fcentres = Vector{Point{dtype}}(nfaces)
    ccentres = Vector{Point{dtype}}(ncells)  
    cvolumes = Vector{dtype}(ncells) 

    # create vector of mesh points
    mesh_points = Point{dtype}[Point(el...) for el in points_data]

    # loop over faces and construct required face data
    mesh_faces = Vector{PolygonalFace{dtype}}(length(faces_data))
    for (faceID, ptsID) in enumerate(faces_data)
        # number of points defining the face
        N = length(ptsID)
        #=
          Select face points for computation of centre and svec. 
          we use tuple because _centre and _svec are methods with 
          multiple implementations, depending on the number of points. 
        =#
        pts = ntuple(i->mesh_points[ptsID[i]], N)
        # check for disaster
        inplane(pts) || error("found non planar face with ID $faceID")
        # now build face        
        f = PolygonalFace{dtype, N}(ptsID,   
                                    _centre(pts),
                                    _svec(pts),
                                    fowners[faceID], 
                                    faceID > length(fneighs) ? 
                                             UInt32(0)      : 
                                             fneighs[faceID])
        mesh_faces[faceID] = f
        # store required data
        fcentres[faceID] = centre(f)
        fsvecs[faceID] = svec(f)
    end

    #=
       Read `owner` and `neighbour` files to obtain a dictionary data 
       structure containing for each cell a list of its faces. 
       Basically, the dictionary `data` will be such that `data[1]` 
       contains the IDs of the faces of cell 1. Reading only the 
       `owner` file is not enough because the dictionary `data`
       would only contain for a given cell the IDs of the faces it 
       owns. By reading the `neighbour` file we also add for each cell
       the IDs of the faces it shares with a neighbour. Note we have 
       Int as the key type and UInt32 as the value.
    =#
    data = DefaultDict(Int, Vector{UInt32}, Vector{UInt32})
    for (faceID, ownerCellID) in enumerate(fowners)
        push!(data[ownerCellID], faceID)
    end
    for (faceID, neighCellID) in enumerate(fneighs)
        push!(data[neighCellID], faceID)
    end

    # Now create `mesh_cells` vector. Note that we need to use indexing 
    # of the `mesh_cells` vector instead of just push! because `data` 
    # is a dictionary which is not ordered. 
    mesh_cells = Vector{PolyHedralCell{dtype}}(length(keys(data)))
    for (cellID, faceIDs) in data
        M = length(faceIDs) # number of faces composing the cell
        fcs = ntuple(i->mesh_faces[faceIDs[i]], M) # tuple of faces
        cell = PolyHedralCell{dtype, M}(tuple(faceIDs...), 
                                        _centreAndVolume(fcs)...)
        mesh_cells[cellID] = cell
        # store required data
        ccentres[cellID] = centre(cell)
        cvolumes[cellID] = volume(cell)
    end

    # create a dict of patches
    patches = (Symbol=>Patch)[k => Patch(k, v...) 
                              for (k, v) in read_boundary(casedir)]
    return Mesh{dtype}(fcentres, fsvecs, fowners, fneighs, 
                       cvolumes, ccentres, patches)
end

end