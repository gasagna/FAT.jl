using Base.Test
using FAT.Meshes

# rotation matrices
Q1 = [1    0          0;
      0    sqrt(3)/2  1/2;
      0    -1/2       sqrt(3)/2]
Q2 = eye(3, 3)      

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

for Q in (Q1, Q2)

    # rotate points
    ab_ = Point(Q*asarray(ab)...)
    bb_ = Point(Q*asarray(bb)...)
    cb_ = Point(Q*asarray(cb)...)
    db_ = Point(Q*asarray(db)...)
    at_ = Point(Q*asarray(at)...)
    bt_ = Point(Q*asarray(bt)...)
    ct_ = Point(Q*asarray(ct)...)
    dt_ = Point(Q*asarray(dt)...)
    nds = [ab_, bb_, cb_, db_, at_, bt_, ct_, dt_]

    # centres and areas of faces
    fcentres = Point{Float64}[]
    fareas = Float64[]

    for faceIDs in faceIDss
        pts = ntuple(i -> nds[faceIDs[i]], 4)
        push!(fcentres, FAT.Meshes._centre(pts...))
        push!(fareas,   FAT.Meshes._area(pts...))
    end

    # build cell
    cc, cv = FAT.Meshes._centreAndVolume(tuple(fareas...), tuple(fcentres...))

    # check centre and volume
    @test distance(cc, Point(Q*asarray(Point(0.5, 0.5, 0.5))...)) < 5e-16
    @test cv ≈ 1.0
end

