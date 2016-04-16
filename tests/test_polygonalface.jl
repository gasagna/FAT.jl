using Base.Test
using FAT.Meshes
# ~~~ Tests for PolygonalFace{Float64, 4} ~~~

# regular face, the unit square
a = Point(0.0, 0.0, 0.0)
b = Point(1.0, 0.0, 0.0)
c = Point(1.0, 1.0, 0.0)
d = Point(0.0, 1.0, 0.0)

# inplane testing
@test inplane((a, b, c)) == true
@test inplane((a, b, c, d)) == true
@test inplane((a, b, c, d + Point(0, 0, 1))) == false

# some tests on triangles
@test FAT.Meshes._centre((a, b, c)) == Point(2, 1, 0)/3
@test FAT.Meshes._centre((c, d, a)) == Point(1, 2, 0)/3
@test FAT.Meshes._svec((a, b, c)) == Point(0.0, 0.0, 1.0)*0.5

# test for quadrilaterals
@test FAT.Meshes._centre((a, b, c, d)) == Point(0.5, 0.5, 0.0)
@test FAT.Meshes._svec((a, b, c, d)) == Point(0.0, 0.0, 1.0)

# test quadface and generic interface
pts = (a, b, c, d)
f = PolygonalFace{Float64, 4}(ntuple(i->UInt32(i), 4), 
             FAT.Meshes._centre(pts), 
             FAT.Meshes._svec(pts), 
             0x00000001,
             0x00000002)

# check area and surface vector
@test area(f) == 1
@test svec(f) == Point(0.0, 0.0, 1.0)
@test ownerID(f) == 0x00000001
@test neighbourID(f) == 0x00000002
@test isInternal(f) == true # neighbourID is not 0
@test isOnBoundary(f) == false # neighbourID is not 0
@test isoutwards(1, f) == true
@test isoutwards(2, f) == false


# Apply some rotation to the points and check area does not vary
# see examples in https://en.wikipedia.org/wiki/Rotation_matrix
Q = [1    0          0;
     0    sqrt(3)/2  1/2;
     0    -1/2       sqrt(3)/2]

ar = Point(Q*asarray(a)...)
br = Point(Q*asarray(b)...)
cr = Point(Q*asarray(c)...)
dr = Point(Q*asarray(d)...)

# test quadface and generic interface
ptsr = (ar, br, cr, dr)
fr = PolygonalFace{Float64, 4}(ntuple(i->UInt32(i), 4), 
              FAT.Meshes._centre(ptsr), 
              FAT.Meshes._svec(ptsr), 
              0x00000001,
              0x00000002)


@test area(fr) ≈ 1
@test svec(fr) == Point(Q*asarray(Point(0.0, 0.0, 1.0))...)

# now multiply points by a constant and move them. 
# Still should be the same.
s = Point(-1, 14, -23)
ptss = (2*ar-s, 2*br-s, 2*cr-s, 2*dr-s)
fs = PolygonalFace{Float64, 4}(ntuple(i->UInt32(i), 4), 
              FAT.Meshes._centre(ptss), 
              FAT.Meshes._svec(ptss), 
              0x00000001,
              0x00000002)


# check area and surface vector
@test area(fs) ≈ 1*2*2
@test distance(svec(fs), Point(Q*asarray(Point(0.0, 0.0, 1.0))...)*2*2) < 1e-14