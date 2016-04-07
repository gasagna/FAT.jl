using Base.Test
using FAT.Meshes

# mainly used for testing and not for actual code
asarray{T}(p::Point{T}) = T[p.x, p.y, p.z]

a = Point(1.0, 0.0, 2.0)
b = Point(0.0, 1.0, 0.0)
c = Point(0.1, 0.2, 0.3)

@test a+b == Point(1.0,  1.0, 2.0)
@test +(a, b, c) == Point(1.1,  1.2, 2.3)
@test -b == Point(0.0,  -1.0, 0.0)
@test a-b == Point(1.0, -1.0, 2.0)
@test a+c == Point(1.1,  0.2, 2.3)
@test a/2 == Point(0.5,  0.0, 1.0)
@test a*b == 0.0
@test a*1 == Point(1.0, 0.0, 2.0)
@test a*0.5 == Point(0.5, 0.0, 1.0)
@test a*0.5 == 0.5*a
@test a*c == 0.7
@test b*c == 0.2

@test norm(a) == sqrt(5)
@test norm(c) == sqrt(0.1^2 + 0.2^2 + 0.3^2)
@test distance(a, b) == sqrt(6)
@test asarray(a) == eltype(a)[1, 0, 2]
@test cross(a, b) + cross(b, a) == Point(0.0, 0.0, 0.0)

# checked with julia cross
@test cross(a, b) == Point(cross(asarray(a), asarray(b))...)
@test cross(a, c) == Point(cross(asarray(a), asarray(c))...)

# check equality in corner case
pminus = Point(0.0, -0.0, 0.0)
pplus  = Point(0.0, +0.0, 0.0)

@test pminus == pplus
@test Point(1, 1, 1) == Point(1, 1, 1)
@test Point(1, 1, -1) != Point(1, 1, 1)
@test Point(1, -1, 1) != Point(1, 1, 1)
@test Point(-1, 1, 1) != Point(1, 1, 1)


