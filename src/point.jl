# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
import Base

export Point,
       distance,
       asarray

"""
    A point in a 3D space. 

    The parameter `T` defines the data type with which
    computations will be performed. This type has three
    fields, with the three spatial coordinates of the point.

"""
struct Point{T<:Real}
    x::T
    y::T
    z::T
end

Point(x::Real, y::Real, z::Real) = Point(promote(x, y, z)...)

Base.show(io::IO, p::Point) = print(io, "($(p.x), $(p.y), $(p.z))")

# mainly used for testing and not for actual code
asarray(p::Point{T}) where {T} = T[p.x, p.y, p.z]

Base.:+(a::Point, b::Point) = Point(a.x+b.x, a.y+b.y, a.z+b.z)
Base.:-(a::Point, b::Point) = Point(a.x-b.x, a.y-b.y, a.z-b.z)
Base.:-(a::Point) = Point(-a.x, -a.y, -a.z)
Base.:/(a::Point, c::Number)  = Point(a.x/c, a.y/c, a.z/c)
Base.:*(a::Point, c::Number)  = Point(a.x*c, a.y*c, a.z*c)
Base.:*(c::Number,  a::Point) = a*c
Base.:*(a::Point, b::Point) = a.x*b.x + a.y*b.y + a.z*b.z
Base.:(==)(a::Point, b::Point) = (a.x==b.x && a.y==b.y && a.z==b.z)

distance(a::Point, b::Point) = norm(a-b)

Base.norm(p::Point) = sqrt(p.x^2 + p.y^2 + p.z^2)
Base.cross(a::Point, b::Point) = Point(a.y*b.z - a.z*b.y, 
                                  a.z*b.x - a.x*b.z, 
                                  a.x*b.y - a.y*b.x)

# make a zero point
Base.zero(::Type{Point{T}}) where {T} = Point(zero(T), zero(T), zero(T))

# eltype
Base.eltype(::Point{T}) where {T} = T