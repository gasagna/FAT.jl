# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
import Base: call

export PolygonalFace, area, centre, svec, ownerID, neighbourID,
       isInternal, isOnBoundary, points_to, isoutwards, inplane,
       QuadFace

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Polygonal face definition ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
immutable PolygonalFace{T<:Real, N}
	points::NTuple{N, UInt32}  # tuple of points ID defining the face. This is better
							   # than a vector, because a tuple is immutable and of fixed size
							   # so its type depends on the input. The parameter N specifies
							   # the number of points defining the face, 3 for triangular, 4 for
							   # quadrilateral, ...
	centre::Point3D{T}         # face centre of mass. The parameter T is used to describe the precision
							   # used for points, centres, and all the rest of the arithmetic.
	svec::Point3D{T}           # surface normal vector. This vector is oriented using the right
							   # hand rule, following OpenFOAM convention. Note that the norm
							   # of this vector is the area of the cell
	owner::UInt32              # ID of the owner cell
	neighbour::UInt32          # ID of the neighbour cell, it is 0 if face is on the boundary.
						       # Note that OpenFoam has -1 but here we have an UInt32 and julia
						       # starts from 1, so 0 is free. Note that we use UInt32 as this would
						       # give a maximum total number of faces equal to 4,294,967,295, which
						       # is already quite big. We could go for UInt64 in a later stage.
end

# Outer constructor. Note that we need to construct the object by providing
# all the necessary information, e.g. centre, ... This is because we do not
# pass the points directly, but only the points ID, and we need to compute
# all necessary information ourselves.
PolygonalFace{T, N}(points::NTuple{N, UInt32},
					centre::Point3D{T},
					svec::Point3D{T},
					owner::UInt32,
					neighbour::UInt32) = PolygonalFace{T, N}(points, centre, svec,
															 owner, neighbour)

# area and face centroid. The area function is not used in real code,
# as we often use the svec function.
@inline area(f::PolygonalFace) = norm(f.svec)
@inline centre(f::PolygonalFace) = f.centre

# get total surface normal vector or each component separately
@inline svec(f::PolygonalFace) = f.svec
@inline svec(f::PolygonalFace, ::Type{Val{:x}}) = f.svec.x
@inline svec(f::PolygonalFace, ::Type{Val{:y}}) = f.svec.y
@inline svec(f::PolygonalFace, ::Type{Val{:z}}) = f.svec.z
@inline svec(f::PolygonalFace, ::Type{Val{1}}) = f.svec.x
@inline svec(f::PolygonalFace, ::Type{Val{2}}) = f.svec.y
@inline svec(f::PolygonalFace, ::Type{Val{3}}) = f.svec.z

# get owner and neighbour cell ID
@inline ownerID(f::PolygonalFace) = f.owner
@inline neighbourID(f::PolygonalFace) = f.neighbour

# is it an internal face or a boundary face
@inline isOnBoundary(f::PolygonalFace) = (f.neighbour == 0x00000000)
@inline isInternal(f::PolygonalFace) = !isOnBoundary(f)

# it seems that for the test case here, the face always points to
# the neighbour and away from its owner, i.e. it points outwards its
# owner. This is not what OpenFOAM Documentation says. What it is said
# is implemented here in code. Hence this could be simplified.
""" Tells whether face points towards `ownerID` or `neighbourID` """
points_to(face::PolygonalFace) = ownerID(face) > neighbourID(face) ? ownerID(face) : neighbourID(face)

""" Tells whether face points outwards the `cellID` cell """
function isoutwards(cellID::Integer, face::PolygonalFace)
	# boundary faces always point outwards
	isOnBoundary(face) && return true
	# if the face points to cellID then clearly is inwards
	return points_to(face) == cellID ? false : true
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for Triangular faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Compute area centroid of triangular face """
@inline _centre{T<:Real}(points::NTuple{3, Point3D{T}}) = +(points...)/3

""" Surface normal vector of a triangular face """
@inline _svec{T<:Real}(points::NTuple{3, Point3D{T}}) = cross(points[2]-points[1],  points[3]-points[2])/2

# TODO speed up this by using an explicit formula
""" Compute area of triangular face as half of the norm of the cross product between any two edges """
@inline _area{T<:Real}(points::NTuple{3, Point3D{T}}) = norm(_svec(points))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for Quadrilateral faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Compute area centroid of quadrilateral face """
@inline function _centre{T<:Real}(points::NTuple{4, Point3D{T}})
	a1 = _area((points[1], points[2], points[3]))
	a2 = _area((points[3], points[4], points[1]))
	c1 = _centre((points[1], points[2], points[3]))
	c2 = _centre((points[3], points[4], points[1]))
	(c1*a1 + c2*a2)/(a1 + a2)
end

""" Compute surface normal vector of quadrilateral face.
	See `Computational Fluid Dynamics: Principles and Applications` Sec. 5.1.2 """
@inline function _svec{T<:Real}(points::NTuple{4, Point3D{T}})
	ΔxA = points[4].x - points[2].x; ΔxB = points[3].x - points[1].x
	ΔyA = points[4].y - points[2].y; ΔyB = points[3].y - points[1].y
	ΔzA = points[4].z - points[2].z; ΔzB = points[3].z - points[1].z
	Point3D{T}((ΔzA*ΔyB - ΔyA*ΔzB)*0.5,
			   (ΔxA*ΔzB - ΔzA*ΔxB)*0.5,
			   (ΔyA*ΔxB - ΔxA*ΔyB)*0.5)
end

""" Area of a quadrilateral """
@inline _area{T<:Real}(points::NTuple{4, Point3D{T}}) = norm(_svec(points))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to check planarity of faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Check all points line in a single plane, by computing the cross
	product of first two neighbouring edges and checking all other edges
	are orthogonal to it """
function inplane{T, N}(points::NTuple{N, Point3D{T}})
	svec = cross(points[2]-points[1], points[3]-points[2])
	for i = 3:N-1
		svec*(points[i+1]-points[i]) ≈ 0.0 || return false
	end
	return true
end

# automatically verified for triangular faces
inplane{T}(points::NTuple{3, Point3D{T}}) = true


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Various aliases used for convenience
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
typealias QuadFace{T} PolygonalFace{T, 4}

# overload call function to allow using QuadFace as a constructor without type parameter
call{T}(q::Type{QuadFace}, points::NTuple{4, UInt32},
					       centre::Point3D{T},
						   svec::Point3D{T},
						   owner::UInt32,
						   neighbour::UInt32) = PolygonalFace{T, 4}(points, centre, svec,
															        owner, neighbour)
