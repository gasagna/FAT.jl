using Base.Test
using FAT.Meshes

casedir = "/Users/davide/Dropbox/Codes/FAT/tests/ldc_test"

m = HexaMesh(casedir)

# we have a cube, with unitary volume and surface area equal to 6
@test_approx_eq_eps sum([volume(cell) for cell in cells(m)])  1 1e-14
@test_approx_eq_eps sum([area(face) for face in filter(isOnBoundary, faces(m))]) 6 1e-10

# sum of all area surface vectors over the boundary patches is zero
@test_approx_eq_eps norm(sum([svec(face) for face in filter(isOnBoundary, faces(m))])) 0 1e-10

# for each cell the sum of the outwards face vectors is zero
fcs = faces(m)
for (cellID, cell) in enumerate(cells(m))
	S = Point3D(0.0, 0.0, 0.0)
	for faceID in facesIDs(cell)
		sign = isoutwards(cellID, fcs[faceID]) ? 1 : -1
		S += svec(fcs[faceID])*sign
	end
	@test norm(S) ≈ 0.0
end

@test npoints(m) == 882
@test nfaces(m) == 1640
@test ncells(m) == 400
@test npatches(m) == 6

@test nboundaryfaces(m) == 880
@test ninternalfaces(m) == nfaces(m) - nboundaryfaces(m)
@test ninternalfaces(m) == 760

@test length(internalfaces(m)) == ninternalfaces(m)
@test length(collect(internalfaces(m))) == ninternalfaces(m)
@test length(boundaryfaces(m)) == nboundaryfaces(m)
@test length(collect(boundaryfaces(m))) == nboundaryfaces(m)

# test enumerate over FaceIterator
for (i, face) in enumerate(internalfaces(m))
	@test isInternal(faces(m)[i]) == true
	@test isOnBoundary(faces(m)[i]) == false
	@test face == faces(m)[i]
end

for (i, face) in enumerate(boundaryfaces(m))
	@test isInternal(faces(m)[i]) == false
	@test isOnBoundary(faces(m)[i]) == true
	@test face == faces(m)[i]
end

for patchname in [:top, :left, :right, :bottom, :back0, :front1]
	for (i, face) in enumerate(faces(m, patchname))
		@test isOnBoundary(face) == true
		@test isInternal(face) == false
		@test face == faces(m)[i]
		p = patch(m, patchname)
		firstfaceID(p) <= i <= lastfaceID(p)
	end
end

for (i, face) in enumerate(faces(m))
	# the first faces are internal, then on the boundary.
	# just look at the lower startface in the patches.
	if i <= ninternalfaces(m)
		@test isInternal(face) == true
		@test isOnBoundary(face) == false
	else
		@test isInternal(face) == false
		@test isOnBoundary(face) == true
	end
end

for face in internalfaces(m)
	# an internal face points from the owner to the neighbour
	# and the ownerID < neighbourID. Unless neighbour is on the boundary.
	# In that case neighbourID = 0
	@test ownerID(face) < neighbourID(face)

	# a internal face always points towards to cell with larger ID 
	@test points_to(face) == neighbourID(face)
end

# for all cells, if a face is on the boundary, it points outwards
for (i, cell) in enumerate(cells(m))
	for faceID in facesIDs(cell)
		face = faces(m)[faceID]
		if isOnBoundary(face)
			@test isoutwards(i, face) == true
		end
	end
end

# test faces on the boundary patches are correctly identified
for (ptchname, coord, value) in zip([:top, :bottom, :left, :right, :back0, :front1], 
									[:y, :y, :x, :x, :z, :z], 
									[1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
	for face in faces(m, ptchname)
		@test getfield(centre(face), coord) == value
	end
end



# test patches
@test patch(m, :top)    == Patch(:top,     20,  761)
@test patch(m, :left)   == Patch(:left,    20,  781)
@test patch(m, :right)  == Patch(:right,   20,  801)
@test patch(m, :bottom) == Patch(:bottom,  20,  821)
@test patch(m, :back0)  == Patch(:back0,  400,  841)
@test patch(m, :front1) == Patch(:front1, 400, 1241)

@test nfaces(patch(m, :top))    == 20
@test nfaces(patch(m, :left))   == 20
@test nfaces(patch(m, :right))  == 20
@test nfaces(patch(m, :bottom)) == 20
@test nfaces(patch(m, :back0))  == 400
@test nfaces(patch(m, :front1)) == 400

@test patchname(patch(m, :top))    == :top
@test patchname(patch(m, :left))   == :left
@test patchname(patch(m, :right))  == :right
@test patchname(patch(m, :bottom)) == :bottom
@test patchname(patch(m, :back0))  == :back0
@test patchname(patch(m, :front1)) == :front1

@test firstfaceID(patch(m, :top))    ==  761
@test firstfaceID(patch(m, :left))   ==  781
@test firstfaceID(patch(m, :right))  ==  801
@test firstfaceID(patch(m, :bottom)) ==  821
@test firstfaceID(patch(m, :back0))  ==  841
@test firstfaceID(patch(m, :front1)) == 1241

@test lastfaceID(patch(m, :top))    ==  780
@test lastfaceID(patch(m, :left))   ==  800
@test lastfaceID(patch(m, :right))  ==  820
@test lastfaceID(patch(m, :bottom)) ==  840
@test lastfaceID(patch(m, :back0))  == 1240
@test lastfaceID(patch(m, :front1)) == 1640


