using Base.Test
using FAT.Meshes

# regular cube
ab = Point3D(0.0, 0.0, 0.0)
bb = Point3D(1.0, 0.0, 0.0)
cb = Point3D(1.0, 1.0, 0.0)
db = Point3D(0.0, 1.0, 0.0)
at = Point3D(0.0, 0.0, 1.0)
bt = Point3D(1.0, 0.0, 1.0)
ct = Point3D(1.0, 1.0, 1.0)
dt = Point3D(0.0, 1.0, 1.0)

nds = [ab, bb, cb, db, at, bt, ct, dt]

# build faces
faceIDss = [(1, 2, 3, 4), # bottom
	        (1, 5, 8, 4), # left
	        (3, 4, 8, 7), # back
	        (2, 6, 7, 3), # right
	        (1, 5, 6, 2), # front
            (5, 6, 7, 8)] # top

fcs = QuadFace{Float64}[]
for faceIDs in faceIDss
	pts = ntuple(i -> nds[faceIDs[i]], 4)
	face = QuadFace{Float64}(faceIDs, 
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
cell = HexaCell(ntuple(i->UInt32(i), 6),
				FAT.Meshes._centreAndVolume(tuple(fcs...))...)

# test volume
@test volume(cell) == 1.0

# test centre
@test centre(cell) == Point3D(0.5, 0.5, 0.5)


# now apply a rotation and check volume is still the same
Q = [1    0          0;
 	 0    sqrt(3)/2  1/2;
 	 0    -1/2       sqrt(3)/2]

ab = Point3D(Q*asarray(ab)...)
bb = Point3D(Q*asarray(bb)...)
cb = Point3D(Q*asarray(cb)...)
db = Point3D(Q*asarray(db)...)
at = Point3D(Q*asarray(at)...)
bt = Point3D(Q*asarray(bt)...)
ct = Point3D(Q*asarray(ct)...)
dt = Point3D(Q*asarray(dt)...)
nds = [ab, bb, cb, db, at, bt, ct, dt]

# do the same
fcs = QuadFace{Float64}[]
for faceIDs in faceIDss
	pts = ntuple(i -> nds[faceIDs[i]], 4)
	face = QuadFace{Float64}(faceIDs, 
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
cell = HexaCell(ntuple(i->UInt32(i), 6),
				FAT.Meshes._centreAndVolume(tuple(fcs...))...)

# test volume
@test volume(cell) == 1.0

# test centre
@test distance(centre(cell), Point3D(Q*asarray(Point3D(0.5, 0.5, 0.5))...)) < 1e-15

# now apply displacement and scaling to the points and verify still works
s = Point3D(-1.0, 2, -4)
nds = [2*ab-s, 2*bb-s, 2*cb-s, 2*db-s, 2*at-s, 2*bt-s, 2*ct-s, 2*dt-s]

# do the same
fcs = QuadFace{Float64}[]
for faceIDs in faceIDss
	pts = ntuple(i -> nds[faceIDs[i]], 4)
	face = QuadFace{Float64}(faceIDs, 
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
cell = HexaCell(ntuple(i->UInt32(i), 6),
				FAT.Meshes._centreAndVolume(tuple(fcs...))...)

# test volume is thrice the scaling
@test volume(cell) ≈ 1.0*2*2*2

# test centre
@test distance(centre(cell), 2*Point3D(Q*asarray(Point3D(0.5, 0.5, 0.5))...) - s) < 1e-15


