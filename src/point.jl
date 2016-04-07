# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
import Base: show,
             norm,
             cross,
             getindex,
             +,
             -,
             /,
             *,
             ==
export Point,
       distance,
       midpoint,
       asarray

"""
    A point in a 3D space. 

    The parameter `T` defines the data type with which
    computations will be performed. This type has three
    fields, with the three spatial coordinates of the point.

"""
immutable Point{T<:Real}
	x::T
	y::T
	z::T
end
Point(x::Real, y::Real, z::Real) = Point(promote(x, y, z)...)
show(io::IO, p::Point) = print(io, "($(p.x), $(p.y), $(p.z))")

# mainly used for testing and not for actual code
asarray{T}(p::Point{T}) = T[p.x, p.y, p.z]

+(a::Point, b::Point) = Point(a.x+b.x, a.y+b.y, a.z+b.z)
-(a::Point, b::Point) = Point(a.x-b.x, a.y-b.y, a.z-b.z)
-(a::Point) = Point(-a.x, -a.y, -a.z)
/(a::Point, c::Number)  = Point(a.x/c, a.y/c, a.z/c)
*(a::Point, c::Number)  = Point(a.x*c, a.y*c, a.z*c)
*(c::Number,  a::Point) = a*c
*(a::Point, b::Point) = a.x*b.x + a.y*b.y + a.z*b.z
==(a::Point, b::Point) = (a.x==b.x && a.y==b.y && a.z==b.z)

norm{T}(p::Point{T}) = sqrt(p.x^2 + p.y^2 + p.z^2)
distance(a::Point, b::Point) = norm(a-b)
cross(a::Point, b::Point) = Point(a.y*b.z - a.z*b.y, 
                                  a.z*b.x - a.x*b.z, 
                                  a.x*b.y - a.y*b.x)
