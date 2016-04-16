using Base.Test
using FAT.Meshes

# regular cube
ab = Point(0.0, 0.0, 0.0)
bb = Point(1.0, 0.0, 0.0)
cb = Point(1.0, 1.0, 0.0)
db = Point(0.0, 1.0, 0.0)
at = Point(0.0, 0.0, 1.0)
bt = Point(1.0, 0.0, 1.0)
ct = Point(1.0, 1.0, 1.0)
dt = Point(0.0, 1.0, 1.0)

nds = [ab, bb, cb, db, at, bt, ct, dt]

# build faces
faceIDss = [(1, 2, 3, 4), # bottom
            (1, 5, 8, 4), # left
            (3, 4, 8, 7), # back
            (2, 6, 7, 3), # right
            (1, 5, 6, 2), # front
            (5, 6, 7, 8)] # top

fcs = PolygonalFace{Float64, 4}[]
for faceIDs in faceIDss
    pts = ntuple(i -> nds[faceIDs[i]], 4)
    face = PolygonalFace{Float64, 4}(faceIDs, 
                             FAT.Meshes._centre(pts),
                             FAT.Meshes._svec(pts),
                             0x00000001, # this is not tested here
                             0x00000002) # this is not tested here
    push!(fcs, face)
end

# test areas are all unitary
for face in fcs
    @test area(face) == 1.0
end

# build cell
cell = PolyHedralCell{Float64, 6}(ntuple(i->UInt32(i), 6),
                FAT.Meshes._centreAndVolume(tuple(fcs...))...)

# test volume
@test volume(cell) == 1.0

# test centre
@test centre(cell) == Point(0.5, 0.5, 0.5)


# now apply a rotation and check volume is still the same
Q = [1    0          0;
     0    sqrt(3)/2  1/2;
     0    -1/2       sqrt(3)/2]

ab = Point(Q*asarray(ab)...)
bb = Point(Q*asarray(bb)...)
cb = Point(Q*asarray(cb)...)
db = Point(Q*asarray(db)...)
at = Point(Q*asarray(at)...)
bt = Point(Q*asarray(bt)...)
ct = Point(Q*asarray(ct)...)
dt = Point(Q*asarray(dt)...)
nds = [ab, bb, cb, db, at, bt, ct, dt]

# do the same
fcs = PolygonalFace{Float64, 4}[]
for faceIDs in faceIDss
    pts = ntuple(i -> nds[faceIDs[i]], 4)
    face = PolygonalFace{Float64, 4}(faceIDs, 
                             FAT.Meshes._centre(pts),
                             FAT.Meshes._svec(pts),
                             0x000000001, 
                             0x000000002)
    push!(fcs, face)
end

# test areas are all unitary
for face in fcs
    @test area(face) ≈ 1.0
end

# build cell
cell = PolyHedralCell{Float64, 6}(ntuple(i->UInt32(i), 6),
                FAT.Meshes._centreAndVolume(tuple(fcs...))...)

# test volume
@test volume(cell) == 1.0

# test centre
@test distance(centre(cell), Point(Q*asarray(Point(0.5, 0.5, 0.5))...)) < 1e-15

# now apply displacement and scaling to the points and verify still works
s = Point(-1.0, 2, -4)
nds = [2*ab-s, 2*bb-s, 2*cb-s, 2*db-s, 2*at-s, 2*bt-s, 2*ct-s, 2*dt-s]

# do the same
fcs = PolygonalFace{Float64, 4}[]
for faceIDs in faceIDss
    pts = ntuple(i -> nds[faceIDs[i]], 4)
    face = PolygonalFace{Float64, 4}(faceIDs, 
                             FAT.Meshes._centre(pts),
                             FAT.Meshes._svec(pts),
                             0x000000001, 
                             0x000000002)
    push!(fcs, face)
end

# test areas are all unitary time twice the scaling
for face in fcs
    @test area(face) ≈ 1.0*2*2
end

# build cell
cell = PolyHedralCell{Float64, 6}(ntuple(i->UInt32(i), 6),
                FAT.Meshes._centreAndVolume(tuple(fcs...))...)

# test volume is thrice the scaling
@test volume(cell) ≈ 1.0*2*2*2

# test centre
@test distance(centre(cell), 2*Point(Q*asarray(Point(0.5, 0.5, 0.5))...) - s) < 1e-15


