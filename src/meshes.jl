# ------------------------------------------------------------------- #
# Copyright 2015-2022, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

import DataStructures: DefaultDict
import HeterogeneousVectors: HVector

import LinearAlgebra

export Patch, 
       patchname,
       firstfaceID,
       lastfaceID,
       nfaces

export Mesh,
       # points functions
       points,
       # patches functions
       patch,
       patches,
       # face functions
       facesvecs,
       facecentres,
       faceownerID,
       faceneighID,
       # nfoo functions
       npatches,
       ncells,
       nfaces,
       npoints,
       # face iterator related
       faceiterator,
       facesIDs,
       nboundaryfaces,
       ninternalfaces,
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
function interpolationWeight(fcentre::Point, ocentre::Point, ncentre::Point, fsvec::Point)
    down = (fcentre - ocentre)*fsvec
    dnei = (ncentre - fcentre)*fsvec 
    return down/(down+dnei)
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
    fsvecs   : face surface vectors
    fcentres : face centres
    fowners  : integer indices representing the owner cell of each face
    fneighs  : integer indices representing the neighbour cell of each face
    cvolumes : cell volumes
    ccentres : cell mass centroid
    patches  : dictionary of patches
    αs : interpolation coefficients for each face
    nboundaryfaces : total number of boundary faces, including 
                     those on empty patches
    ninternalfaces : the total number of internal faces
"""
struct Mesh{T}
    points::Vector{Point{T}}
    fsvecs::Vector{Point{T}}
    fcentres::Vector{Point{T}}
    fowners::Vector{UInt32}
    fneighs::Vector{UInt32}
    cvolumes::Vector{T}
    ccentres::Vector{Point{T}}
    patches::Dict{Symbol, Patch}
    αs::Vector{T}
    nboundaryfaces::Int
    ninternalfaces::Int
    function Mesh{T}(points, fsvecs, fcentres, fowners, fneighs, 
                     cvolumes, ccentres, patches) where {T}
        # number of boundary faces
        nboundaryfaces = sum([nfaces(v) for (p, v) in patches])
        # the remaining are internal faces, for which we need an α
        ninternalfaces = length(fcentres) - nboundaryfaces
        # calculate interpolation coefficients
        αs = Vector{T}(undef, ninternalfaces)
        for i = 1:ninternalfaces
            αs[i] = interpolationWeight(fcentres[i], 
                                        ccentres[fowners[i]], 
                                        ccentres[fneighs[i]], fsvecs[i])
        end
        new{T}(points, fsvecs, fcentres, fowners, fneighs, cvolumes, 
            ccentres, patches, αs, nboundaryfaces, ninternalfaces)
    end
end

# data type of the mesh points, areas, volumes, ...
Base.eltype(::Mesh{T}) where {T} = T

function Base.show(io::IO, mesh::Mesh; gap::AbstractString=" ")
    print(io, "Mesh object \n")
    print(io, "$gap ~ $(ncells(mesh)) cells                 \n")
    print(io, "$gap ~ $(nfaces(mesh)) total faces           \n")
    print(io, "$gap ~ $(ninternalfaces(mesh)) internal faces\n")
    print(io, "$gap ~ $(nboundaryfaces(mesh)) boundary faces\n")
    print(io, "$gap ~ $(length(keys(patches(mesh)))) patches: \n")
    for (patchname, v) in mesh.patches
        print(io, "$gap   ~ $(string(patchname)): $(nfaces(v)) faces\n")
    end
end

# include face iterator
include("faceiterator.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~
# === Getter functions ===
# ~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ cell functions ~~~

" Total number of cells in the mesh "
ncells(m::Mesh) = length(m.cvolumes)

" Get cell volumes "
cellvolumes(m::Mesh) = m.cvolumes

" Get `dir` coordinates of cell centres "
function cellcentres(m::Mesh{T}, dir::Symbol) where {T}
    out = Vector{T}(undef, ncells(m))
    cc = m.ccentres
    for i in 1:length(out)
        out[i] = getfield(cc[i], dir)
    end
    out
end

" Get cell centres "
cellcentres(m::Mesh) = 
    (cellcentres(m, :x), cellcentres(m, :y), cellcentres(m, :z))

# ~~~ face functions ~~~

" Total number of faces in the mesh "
nfaces(m::Mesh) = length(m.fsvecs)

" Face surface vectors "
facesvecs(m::Mesh) = m.fsvecs

" Face centres "
facecentres(m::Mesh) = m.fcentres

" Face cell owner id "
faceownerID(m::Mesh) = m.fowners

" Face cell neighbour id "
faceneighID(m::Mesh) = m.fneighs

" Number of faces on the domain boundary "
@inline nboundaryfaces(m::Mesh) = m.nboundaryfaces

" Number of internal faces "
@inline ninternalfaces(m::Mesh) = m.ninternalfaces

# ~~~ patch functions ~~~

" Total number of patches in the mesh "
npatches(m::Mesh) = length(m.patches)

" Dictionary of mesh patches "
patches(m::Mesh) = m.patches

" Get a patch named `pname`"
patch(m::Mesh, pname::Symbol) = m.patches[pname]

# ~~~ points functions ~~~

" Total number of points in the mesh "
npoints(m::Mesh) = length(m.points)

" Points in the mesh "
points(m::Mesh) = m.points

# ~~~~~~~~~~~~~~~~~~~~~~~~
# === Mesh reader ===
# ~~~~~~~~~~~~~~~~~~~~~~~~

""" 
        Mesh(casedir, mtype)

    Read the mesh in an OpenFoam case directory 

    Parameters
    ----------
    casedir : the OpenFOAM case directory
    mtype   : [optional - default Float64] datatype for points, centres, 
              volumes, ...

"""
function Mesh(casedir::AbstractString, ::Type{T}=Float64, tol::Real=1e-15) where {T<:Real}
    # check before loading
    iscasedir(casedir) || error("$casedir is not an OpenFoam case!")

    # create vector of mesh points. It he
    points_data = reader(casedir, "points", T)
    points = Point{T}[Point(el...) for el in points_data]
    
    # read face information
    faces_data = reader(casedir, "faces")
    local nfaces = length(faces_data)
    
    # we also need to construct these two along with the faces
    fcentres = Vector{Point{T}}(undef, nfaces)
    fsvecs   = Vector{Point{T}}(undef, nfaces)
    fareas   = Vector{T}(undef, nfaces)

    # Temporary vector for storing points in a faces
    pts = Vector{Point{T}}(undef, 10)

    # Now build face properties 
    for faceID = 1:nfaces
        # IDs of points making the current face
        ptsIDs = faces_data[faceID]
        # copy points making current face into temporary pts
        for (i, ptsID) in enumerate(ptsIDs)
            pts[i] = points[ptsID]
        end
        # check for disaster
        _inplane(pts, length(ptsIDs))#, tol) || error("found non planar face with ID $faceID")
        fcentres[faceID] = _centre(pts, length(ptsIDs))
        fsvecs[faceID]   = _svec(pts, length(ptsIDs))
        fareas[faceID]   = LinearAlgebra.norm(fsvecs[faceID])
    end

    # read cell owner information
    fowners = reader(casedir, "owner")::Vector{UInt32}
    fneighs = reader(casedir, "neighbour")::Vector{UInt32}

    # size of mesh
    local ncells = max(maximum(fowners), maximum(fneighs))

    # will construct this as well
    ccentres = Vector{Point{T}}(undef, ncells)  
    cvolumes = Vector{T}(undef, ncells) 

    #=
       Read `owner` and `neighbour` files to obtain a vector of vectors, 
       in which the `i`-th vector contains a list of faces forming cell `i`.
       Reading only the `owner` file is not enough because the resulting
       data structure would only contain, for a given cell, the IDs of the 
       faces it owns. By reading the `neighbour` file we also add, for each 
       cell, the IDs of the faces it shares with a neighbour.
    =#
    data = Vector{UInt32}[UInt32[] for i = 1:ncells]
    for fcellIDs in (fowners, fneighs)
        for faceID in eachindex(fcellIDs)
            push!(data[fcellIDs[faceID]], faceID)
        end
    end

    # Temporaries for constructing the cell information
    areas = Vector{T}(undef, 10)
    centres = Vector{Point{T}}(undef, 10)

    # Now construct cell information
    for (cellID, faceIDs) in enumerate(data)  # for each cell
        for (i, faceID) in enumerate(faceIDs) # for each face in the cell
            areas[i] = fareas[faceID]         # collect face area
            centres[i] = fcentres[faceID]     # collect face centre
        end
        # then perform computation
        ccentres[cellID], cvolumes[cellID] = (
            _centreAndVolume(areas, centres, length(faceIDs)) )
    end

    # create a dict of patches
    patches = Dict{Symbol, Patch}(k => Patch(k, v...) 
                              for (k, v) in read_boundary(casedir))
    return Mesh{T}(points, fsvecs, fcentres, fowners, fneighs, 
                       cvolumes, ccentres, patches)
end