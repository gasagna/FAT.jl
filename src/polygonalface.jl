# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for Triangular faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Compute area centroid of triangular face """
@inline _centre{T<:Real}(points::NTuple{3, Point{T}}) = +(points...)/3

""" Surface normal vector of a triangular face """
@inline _svec{T<:Real}(points::NTuple{3, Point{T}}) = 
    cross(points[2]-points[1],  points[3]-points[2])/2

# TODO speed up this by using an explicit formula
""" Compute area of triangular face as half of the norm of the 
    cross product between any two edges """
@inline _area{T<:Real}(points::NTuple{3, Point{T}}) = 
    norm(_svec(points))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for Quadrilateral faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Compute area centroid of quadrilateral face """
@inline function _centre{T<:Real}(points::NTuple{4, Point{T}})
    a1 = _area((points[1], points[2], points[3]))
    a2 = _area((points[3], points[4], points[1]))
    c1 = _centre((points[1], points[2], points[3]))
    c2 = _centre((points[3], points[4], points[1]))
    (c1*a1 + c2*a2)/(a1 + a2)
end

""" Compute surface normal vector of quadrilateral face.
    See `Computational Fluid Dynamics: Principles and 
    Applications` Sec. 5.1.2 
"""
@inline function _svec{T<:Real}(points::NTuple{4, Point{T}})
    ΔxA = points[4].x - points[2].x; ΔxB = points[3].x - points[1].x
    ΔyA = points[4].y - points[2].y; ΔyB = points[3].y - points[1].y
    ΔzA = points[4].z - points[2].z; ΔzB = points[3].z - points[1].z
    Point{T}((ΔzA*ΔyB - ΔyA*ΔzB)*0.5,
               (ΔxA*ΔzB - ΔzA*ΔxB)*0.5,
               (ΔyA*ΔxB - ΔxA*ΔyB)*0.5)
end

""" Area of a quadrilateral """
@inline _area{T<:Real}(points::NTuple{4, Point{T}}) = 
    norm(_svec(points))

# need to implement method for generic number of faces

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to check planarity of faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Check all points line in a single plane, by computing the cross
    product of first two neighbouring edges and checking all other 
    edges are orthogonal to it 
"""
function _inplane{T, N}(points::NTuple{N, Point{T}})
    svec = cross(points[2]-points[1], points[3]-points[2])
    for i = 3:N-1
        svec*(points[i+1]-points[i]) ≈ 0.0 || return false
    end
    return true
end

# automatically verified for triangular faces
_inplane{T}(points::NTuple{3, Point{T}}) = true