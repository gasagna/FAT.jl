# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

for fname in [:_centre, :_svec, :_inplane]
    @eval function $fname{T<:Real}(pts::Vector{Point{T}}, Npts::Integer)
              Npts == 3 && return $fname(pts[1], pts[2], pts[3])
              Npts == 4 && return $fname(pts[1], pts[2], pts[3], pts[4])
              error("do not know how to compute $fname of such a face")
          end
end 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for Triangular faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Compute area centroid of triangular face """
@inline _centre{T<:Real}(pt1::Point{T}, pt2::Point{T}, pt3::Point{T}) =
     (pt1 + pt2 + pt3)/3

""" Surface normal vector of a triangular face """
@inline _svec{T<:Real}(pt1::Point{T}, pt2::Point{T}, pt3::Point{T}) =
    cross(pt2-pt1, pt3-pt2)/2

# TODO speed up this by using an explicit formula
""" Compute area of triangular face as half of the norm of the 
    cross product between any two edges """
@inline _area{T<:Real}(pt1::Point{T}, pt2::Point{T}, pt3::Point{T}) =
    norm(_svec(pt1, pt2, pt3))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for Quadrilateral faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Compute area centroid of quadrilateral face """
@inline function _centre{T<:Real}(pt1::Point{T}, 
                                  pt2::Point{T}, 
                                  pt3::Point{T}, 
                                  pt4::Point{T})
    a1 =   _area(pt1, pt2, pt3)
    a2 =   _area(pt3, pt4, pt1)
    c1 = _centre(pt1, pt2, pt3)
    c2 = _centre(pt3, pt4, pt1)
    (c1*a1 + c2*a2)/(a1 + a2)
end

""" Compute surface normal vector of quadrilateral face.
    See `Computational Fluid Dynamics: Principles and 
    Applications` Sec. 5.1.2 
"""
@inline function _svec{T<:Real}(pt1::Point{T}, 
                                pt2::Point{T}, 
                                pt3::Point{T}, 
                                pt4::Point{T})
    ΔxA = pt4.x - pt2.x; ΔxB = pt3.x - pt1.x
    ΔyA = pt4.y - pt2.y; ΔyB = pt3.y - pt1.y
    ΔzA = pt4.z - pt2.z; ΔzB = pt3.z - pt1.z
    Point{T}((ΔzA*ΔyB - ΔyA*ΔzB)*0.5,
               (ΔxA*ΔzB - ΔzA*ΔxB)*0.5,
               (ΔyA*ΔxB - ΔxA*ΔyB)*0.5)
end

""" Area of a quadrilateral """
@inline _area{T<:Real}(pt1::Point{T}, 
                       pt2::Point{T}, 
                       pt3::Point{T}, 
                       pt4::Point{T}) = norm(_svec(pt1, pt2, pt3, pt4))

# need to implement method for generic number of faces

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to check planarity of faces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" 
    Check all pts line in a single plane, by computing the cross
    product of first two neighbouring edges and checking all other 
    edges are orthogonal to it.
"""
function _inplane{T}(pt1::Point{T}, 
                     pt2::Point{T}, 
                     pt3::Point{T}, 
                     pt4::Point{T})
    return cross(pt2-pt1, pt3-pt2)*(pt4-pt3) ≈ 0.0 ? true : false
end

_inplane{T}(pt1::Point{T}, pt2::Point{T}, pt3::Point{T}) = true