using Base.Test
using FAT.Meshes

# four points on the unit square 
a = Point(0.0, 0.0, 0.0)
b = Point(1.0, 0.0, 0.0)
c = Point(1.0, 1.0, 0.0)
d = Point(0.0, 1.0, 0.0)

# _inplane testing
@test FAT.Meshes._inplane(a, b, c) == true
@test FAT.Meshes._inplane(a, b, c, d) == true
@test FAT.Meshes._inplane(a, b, c, d + Point(0, 0, 1)) == false

# ~~~ tests for triangles ~~~
@test FAT.Meshes._centre(a, b, c) == Point(2, 1, 0)/3
@test FAT.Meshes._centre(c, d, a) == Point(1, 2, 0)/3

# the surface vector is oriented using the right hand rule.
@test FAT.Meshes._svec(a, b, c) == Point(0.0, 0.0, 1.0)*0.5
@test FAT.Meshes._svec(c, d, a) == Point(0.0, 0.0, 1.0)*0.5
@test FAT.Meshes._svec(a, d, c) == -FAT.Meshes._svec(c, d, a)
@test FAT.Meshes._svec(a, d, c) ==  FAT.Meshes._svec(d, c, a)
@test FAT.Meshes._svec(d, c, a) ==  FAT.Meshes._svec(c, a, d)

# ~~~ test for quadrilaterals ~~~
@test FAT.Meshes._centre(a, b, c, d) == Point(0.5, 0.5, 0.0)

@test FAT.Meshes._svec(a, b, c, d) == Point(0.0, 0.0, 1.0)
@test FAT.Meshes._svec(b, c, d, a) == Point(0.0, 0.0, 1.0)
@test FAT.Meshes._svec(c, d, a, b) == Point(0.0, 0.0, 1.0)
@test FAT.Meshes._svec(d, a, b, c) == Point(0.0, 0.0, 1.0)
@test FAT.Meshes._svec(c, b, a, d) == - FAT.Meshes._svec(d, a, b, c)