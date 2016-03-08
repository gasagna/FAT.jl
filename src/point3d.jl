# -------------------------------------------------------------- #
# Copyright 2015, Davide Lasagna, AFM, University of Southampton #
# -------------------------------------------------------------- #
import Base: show, norm, cross, +, -, /, *, ==
export Point3D, distance, midpoint, asarray


immutable Point3D{T<:Real}
	x::T
	y::T
	z::T
end
Point3D(x::Real, y::Real, z::Real) = Point3D(promote(x, y, z)...)
show(io::IO, p::Point3D) = print(io, "($(p.x), $(p.y), $(p.z))")

# mainly used for testing and not for actual code
asarray{T}(p::Point3D{T}) = T[p.x, p.y, p.z]

+(a::Point3D, b::Point3D) = Point3D(a.x+b.x, a.y+b.y, a.z+b.z)
-(a::Point3D, b::Point3D) = Point3D(a.x-b.x, a.y-b.y, a.z-b.z)
-(a::Point3D) = Point3D(-a.x, -a.y, -a.z)
/(a::Point3D, c::Number)  = Point3D(a.x/c, a.y/c, a.z/c)
*(a::Point3D, c::Number)  = Point3D(a.x*c, a.y*c, a.z*c)
*(c::Number,  a::Point3D) = a*c
# dot product
*(a::Point3D, b::Point3D) = a.x*b.x + a.y*b.y + a.z*b.z
==(a::Point3D, b::Point3D) = (a.x==b.x && a.y==b.y && a.z==b.z)

norm{T}(p::Point3D{T}) = sqrt(p.x^2 + p.y^2 + p.z^2)
midpoint(a::Point3D, b::Point3D) = (a+b)/2
distance(a::Point3D, b::Point3D) = norm(a-b)
cross(a::Point3D, b::Point3D) = Point3D(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)
