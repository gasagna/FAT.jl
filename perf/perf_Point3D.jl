import Base: +
using FAT.Meshes

const a = Point3D(1.0, 0.0, 2.0)
const b = Point3D(0.0, 1.0, 0.0)
const c = Point3D(0.1, 0.2, 0.3)
const points = (a, b, c)


# +(a::Point3D, b::Point3D, c::Point3D) = Point3D(a.x+b.x+c.x, a.y+b.y+c.y, a.z+b.z+c.z)

@inline test{T<:Real}(points::NTuple{3, Point3D{T}}) = +(points...)/3


function t1(points)
	p = points[1]
	for i = 1:100000000
		p = p + test(points)
	end
	p
end

@time t1(points)

println(minimum([@elapsed t1(points) for i = 1:100]))