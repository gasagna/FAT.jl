# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Meshes

import DataStructures: DefaultDict
import Base: start, next, done, length, eltype, show
import FAT.OFIO: reader, read_boundary, isCaseDir


export Patch, patchname, firstfaceID, lastfaceID, nfaces
export Mesh, HexaMesh, cells, faces, patch, points, npatches, patches,
	   ncells, nfaces, npoints, nboundaryfaces, ninternalfaces,
	   internalfaces, boundaryfaces, cellcentres, enumerateFaces

include("point3d.jl")
include("polygonalface.jl")
include("polyhedralcell.jl")
include("patches.jl")

""" Calculate interpolation weights from cell centres to face centre 

	Notes
	-----
	See formula 6.31 in the book The Finite Volume Method in Computational 
	Fluid Dynamics: An Advanced Introduction with OpenFOAM® and Matlab. See 
	figure 6.24
"""
function interpolationWeight(owner::PolyHedralCell, neighbour::PolyHedralCell, face::PolygonalFace)
    down = (centre(face) - centre(owner))*svec(face)
    dnei = (centre(neighbour) - centre(face))*svec(face)
    down/(down+dnei)
end

# ~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Mesh definition ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~

abstract Mesh{T<:Real} 

# data type of the mesh points, faces, cells, ...
eltype{T}(::Mesh{T}) = T


""" Type for meshes with hexahedral cells 

	The type parameters are the type of the data, i.e. points, and 
	the type of the integer used to represent the neighbouring cell
	ID in the cells vector.

	Fields
	------
	nboundaryfaces : the total number of boundary faces, including those on empty patches
	ninternalfaces : the total number of internal faces
"""
immutable HexaMesh{T} <: Mesh{T}
	points::Vector{Point3D{T}}
	faces::Vector{QuadFace{T}}
	patches::Dict{Symbol, Patch}
	cells::Vector{PolyHedralCell{T, 6}}
	αs::Vector{T}
	nboundaryfaces::Int
	ninternalfaces::Int
	function HexaMesh(points, faces, patches, cells)
		# 
		nboundaryfaces = sum([nfaces(v) for (p, v) in patches])
		ninternalfaces = length(faces) - nboundaryfaces
		# Pre compute the interpolation weights for each internal face. This is a parameter used for linear interpolation 
		# of the internal face values from the cell centres of the owner and neighbour. Another solution would
		# be to have the faces have an alpha field and use that. However, the problem is that the faces are 
		# immutable types and you do not know the cell centres until you have created all the faces. Using a 
		# mutable type for faces significantly reduce performance.
		αs = Vector{T}(ninternalfaces)
		for i = 1:ninternalfaces
			αs[i] = interpolationWeight(cells[ownerID(faces[i])], cells[neighbourID(faces[i])], faces[i])
		end
		new(points, faces, patches, cells, αs, nboundaryfaces, ninternalfaces)
	end
end

# Generic iterator over faces
immutable FaceIterator{I}
	facesVec::I
	first::Int # index of the first
	last::Int  # index of the last
	function FaceIterator(facesVec::I, first::Int, last::Int)
		last >= first && first >=1 && last <= length(facesVec) || error("wrong face iterator specification")
		new(facesVec, first, last)
	end
end
# efficient iterators over the internal and boundary faces
internalfaces(m::HexaMesh) = FaceIterator{typeof(m.faces)}(m.faces, 1, m.ninternalfaces)
boundaryfaces(m::HexaMesh) = FaceIterator{typeof(m.faces)}(m.faces, m.ninternalfaces+1, nfaces(m))


# iterate over the faces of a particular patch
faces(m::HexaMesh, p::Symbol) = 
	FaceIterator{typeof(m.faces)}(m.faces, 
								  firstfaceID(m.patches[p]), 
								  lastfaceID(m.patches[p]))

start(fi::FaceIterator) = fi.first
next(fi::FaceIterator, i::Int) = fi.facesVec[i], i+1
done(fi::FaceIterator, i::Int) = i == fi.last + 1
length(fi::FaceIterator) = fi.last - fi.first + 1
eltype{I}(::Type{FaceIterator{I}}) = eltype(I)

# Enumerate over the faces of the FaceIterator.
# This will return the appropriate face ID and the face, as using the 
# standard enumerate would start from 1. Do not get confused.
# the following lines are copied from iterator.jl in Julia Base
immutable EnumerateFaces{I}
    faceitr::I
end
enumerateFaces(faceitr::FaceIterator) = EnumerateFaces(faceitr)

length(e::EnumerateFaces) = length(e.faceitr)
start(e::EnumerateFaces) = (e.faceitr.first, start(e.faceitr))
function next(e::EnumerateFaces, state)
    n = next(e.faceitr, state[2])
    (state[1], n[1]), (state[1]+1, n[2])
end
done(e::EnumerateFaces, state) = done(e.faceitr, state[2])
eltype{I}(::Type{EnumerateFaces{I}}) = Tuple{Int, eltype(I)}

# This one uses directly the faces vector for efficiency
faces(m::HexaMesh) = m.faces
nfaces(m::HexaMesh) = length(m.faces)

# number of boundary and internal faces 
nboundaryfaces(m::HexaMesh) = m.nboundaryfaces
ninternalfaces(m::HexaMesh) = m.ninternalfaces

# number of mesh points and points getter
npoints(m::HexaMesh) = length(m.points)
points(m::HexaMesh) = m.points

# number of cells and cells getter
ncells(m::HexaMesh) = length(m.cells)
cells(m::HexaMesh) = m.cells

# number of patches and patch getter
npatches(m::HexaMesh) = length(m.patches)
patches(m::HexaMesh) = m.patches
patch(m::HexaMesh, p::Symbol) = m.patches[p]

# get cell centres
function cellcentres(m::HexaMesh, dir::Symbol)
	out = Vector{Float64}(ncells(m))
	for (i, cell) in enumerate(cells(m))
		out[i] = getfield(centre(cell), dir)
	end
	out
end

cellcentres(m::HexaMesh) = 
	(cellcentres(m, :x), cellcentres(m, :y), cellcentres(m, :z))

function show(io::IO, mesh::HexaMesh; space::AbstractString=" ")
	print(io, "Hexahedral mesh object at $(object_id(mesh)):  \n")
	print(io, "$space ~ $(ncells(mesh)) cells                 \n")
	print(io, "$space ~ $(npoints(mesh)) points               \n")
	print(io, "$space ~ $(nfaces(mesh)) total faces           \n")
	print(io, "$space ~ $(nboundaryfaces(mesh)) boundary faces\n")
	print(io, "$space ~ $(npatches(mesh)) patches: \n")
	for (patchname, v) in mesh.patches
		print(io, "$space   ~ $(string(patchname)): $(nfaces(v)) faces\n")
	end
end


""" Read the mesh in an OpenFoam case directory 

	Parameters
	==========
	casedir  : the case directory
	RealType : the type of the space and flow variables. Defaults to Float64

	Returns
	=======
	HexaMesh : the mesh 
"""
function HexaMesh{Tf<:Real}(casedir::AbstractString, RealType::Type{Tf}=Float64)

	# check before loading
	isCaseDir(casedir) || error("$casedir is not an OpenFoam case directory!")

	# read all data. These are 1-based, and not 0-based as the OF files
	points_data = reader(casedir, "points")
	faces_data  = reader(casedir, "faces")
	ownerOf     = reader(casedir, "owner")
	neighbourOf = reader(casedir, "neighbour")

	# create vector of Point3D
	mesh_points = Point3D{RealType}[Point3D(el...) for el in points_data]

	# create vector of QuadFaces
	mesh_faces = Vector{QuadFace{RealType}}(length(faces_data))
	for (faceID, ptsID) in enumerate(faces_data)
		# select face points for computation of centre and svec. 
		# we use tuple because _centre and _svec are methods with multiple
		# implementations, depending on the number of points. We might want
		# to adopt a single general formula that works for any number of points
		# avoiding the need to create a tuple of Point3D every time.
		# check all points are on the same plane
		pts = ntuple(i->mesh_points[ptsID[i]], 4)
		inplane(pts) || error("found non planar face with ID $faceID")
		# now build face		
		f = QuadFace{RealType}(ptsID,   
			 		  		   _centre(pts),
					  		   _svec(pts),
					  		   ownerOf[faceID], 
					  		   faceID > length(neighbourOf) ? 
					  		  		    UInt32(0)           : 
					  		  		    neighbourOf[faceID])
		mesh_faces[faceID] = f
	end

	# Read `owner` and `neighbour` files to obtain a dictionary data structure
	# containing for each cell a list of its faces. Basically, the dictionary 
	# `data` will be such that `data[1]` will contain the IDs of the faces of cell 1.
	# Reading only the `owner` file is not enough because the dictionary `data`
	# would only contain for a given cell the IDs of the faces it owns. By reading
	# the `neighbour` file we also add for each cell the IDs of the faces it shares
	# with a neighbour.
	# Note we have Int as the key type and UInt32 as the value
	data = DefaultDict(Int, Vector{UInt32}, Vector{UInt32})
	for (faceID, ownerCellID) in enumerate(ownerOf)
		push!(data[ownerCellID], faceID)
	end
	for (faceID, neighCellID) in enumerate(neighbourOf)
		push!(data[neighCellID], faceID)
	end

	# Now create `mesh_cells` vector. Note that we need to use indexing of the 
	# `mesh_cells` vector instead of just push! because `data` is a dictionary
	# which is not ordered. 
	mesh_cells = Vector{PolyHedralCell{RealType, 6}}(length(keys(data)))
	for (cellID, faceIDs) in data
		fcs = ntuple(i->mesh_faces[faceIDs[i]], 6) 
		mesh_cells[cellID] = PolyHedralCell{RealType, 6}(tuple(faceIDs...), _centreAndVolume(fcs)...)
	end

	# create a dict of patches
	patches = (Symbol=>Patch)[k => Patch(k, v...) for (k, v) in read_boundary(casedir)]
	return HexaMesh{RealType}(mesh_points, mesh_faces, patches, mesh_cells)
end

end

























