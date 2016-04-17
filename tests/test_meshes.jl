using Base.Test
using FAT.Meshes

casedir = "./ldc_test"
m = Mesh(casedir)

# ~~~~~~~~~~~~~~~~~~~
# TEST CELL FUNCTIONS
# ~~~~~~~~~~~~~~~~~~~
@test ncells(m) == 400

# we have an unitary cube
@test abs(sum(m.cvolumes) - 1) < 1e-14

# all cells should have equal volume. We round here because of finite 
# precision arithmetic when loading the data files and in the computations
@test all(round(m.cvolumes, 6) .== 1.0/400)

# test owner ID is always lower then neighbour ID
fo = m.fowners
fn = m.fneighs
for faceID in 1:ninternalfaces(m)
    @test fo[faceID] < fn[faceID]
end

# cell centres is not tested as only used for plotting

# for each cell the sum of the outwards face vectors is zero
# we currently do not hold the cell->faces information, but 
# only the other way round. Hence loop over faces here:
fs = facesvecs(m)
fo = m.fowners
fn = m.fneighs
s = zeros(Point{Float64}, ncells(m))
for faceID in 1:nfaces(m)
    # if it is an internal face we need to consider the neighbour cell
    if faceID <= ninternalfaces(m)
        s[fn[faceID]] -= fs[faceID]
    end
    s[fo[faceID]] += fs[faceID]
end
@test maximum(map(norm, s)) == 0.0

# ~~~~~~~~~~~~~~~~~~~
# TEST FACE FUNCTIONS
# ~~~~~~~~~~~~~~~~~~~
@test nfaces(m) == 1640
@test nboundaryfaces(m) == 880
@test ninternalfaces(m) == nfaces(m) - nboundaryfaces(m)
@test ninternalfaces(m) == 760

# sum of all area surface vectors over the boundary patches is zero,
# while total area is six
fs = facesvecs(m)
S = Point(0.0, 0.0, 0.0)
A = 0.0
for (patchname, patch) in patches(m) 
    S += sum(fs[facesIDs(patch)])
    A += sum(map(norm, fs[facesIDs(patch)]))
end
@test S == Point(0.0, 0.0, 0.0)
@test A ≈ 6

# test faces on the boundary patches are correctly 
fc = facecentres(m)
for (ptch, coord, value) in zip([:top, :bottom, :left, :right, :back0, :front1], 
                                    [:y, :y, :x, :x, :z, :z], 
                                    [1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
      for faceID in facesIDs(m, ptch)
          @test getfield(fc[faceID], coord) == value
      end
end

@test length(facesIDs(m, :internal)) == ninternalfaces(m)
@test length(facesIDs(m, :boundary)) == nboundaryfaces(m)
@test ninternalfaces(m) + nboundaryfaces(m) == nfaces(m)

# ~~~~~~~~~~~~~~~~~~~~
# TEST PATCH FUNCTIONS
# ~~~~~~~~~~~~~~~~~~~~
@test npatches(m) == 6
@test sort(collect(keys(patches(m)))) == sort([:top, :left, :right, 
                                               :bottom, :back0, :front1])
@test patch(m, :top)    == Patch(:top,    false,  20,  761)
@test patch(m, :left)   == Patch(:left,   false,  20,  781)
@test patch(m, :right)  == Patch(:right,  false,  20,  801)
@test patch(m, :bottom) == Patch(:bottom, false,  20,  821)
@test patch(m, :back0)  == Patch(:back0,  true,  400,  841)
@test patch(m, :front1) == Patch(:front1, true,  400, 1241)

@test isempty(patch(m, :top))    == false
@test isempty(patch(m, :left))   == false
@test isempty(patch(m, :right))  == false
@test isempty(patch(m, :bottom)) == false
@test isempty(patch(m, :back0))  == true
@test isempty(patch(m, :front1)) == true

@test patchname(patch(m, :top))    == :top
@test patchname(patch(m, :left))   == :left
@test patchname(patch(m, :right))  == :right
@test patchname(patch(m, :bottom)) == :bottom
@test patchname(patch(m, :back0))  == :back0
@test patchname(patch(m, :front1)) == :front1

@test nfaces(patch(m, :top))    == 20
@test nfaces(patch(m, :left))   == 20
@test nfaces(patch(m, :right))  == 20
@test nfaces(patch(m, :bottom)) == 20
@test nfaces(patch(m, :back0))  == 400
@test nfaces(patch(m, :front1)) == 400

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

@test facesIDs(patch(m, :top))    ==  761:780
@test facesIDs(patch(m, :left))   ==  781:800
@test facesIDs(patch(m, :right))  ==  801:820
@test facesIDs(patch(m, :bottom)) ==  821:840
@test facesIDs(patch(m, :back0))  ==  841:1240
@test facesIDs(patch(m, :front1)) == 1241:1640

# ~~~~~~~~~~~~~~~~~~
# TEST FACE ITERATOR 
# ~~~~~~~~~~~~~~~~~~

# test facesIDs
@test facesIDs(m, :internal) ==    1:760
@test facesIDs(m, :boundary) ==  761:1640
@test facesIDs(m, :top)      ==  761:780
@test facesIDs(m, :left)     ==  781:800
@test facesIDs(m, :right)    ==  801:820
@test facesIDs(m, :bottom)   ==  821:840
@test facesIDs(m, :back0)    ==  841:1240
@test facesIDs(m, :front1)   == 1241:1640

# test faceiterator
@test faceiterator(m, :internal) == zip(   1:760  ,    1:760         )
@test faceiterator(m, :boundary) == zip( 761:1640 , ( 761:1640) - 760)
@test faceiterator(m, :top)      == zip( 761:780  , ( 761:780 ) - 760)
@test faceiterator(m, :left)     == zip( 781:800  , ( 781:800 ) - 760)
@test faceiterator(m, :right)    == zip( 801:820  , ( 801:820 ) - 760)
@test faceiterator(m, :bottom)   == zip( 821:840  , ( 821:840 ) - 760)
@test faceiterator(m, :back0)    == zip( 841:1240 , ( 841:1240) - 760)
@test faceiterator(m, :front1)   == zip(1241:1640 , (1241:1640) - 760)

for faceID in facesIDs(m, :internal)
  # an internal face points from the owner to the neighbour and so 
  # ownerID < neighbourID. Only internal faces have a neighbour. 
  @test m.fowners[faceID] < m.fneighs[faceID]
end

# ~~~~~~~~~~~~~~~~~~~~
# TEST INTERNAL FIELDS 
# ~~~~~~~~~~~~~~~~~~~~
#  ~~~ points

@test npoints(m) == 882

# they are in the same order as in the file
@test points(m)[1]     == Point(0.00, 0, 0)
@test points(m)[2]     == Point(0.05, 0, 0)
@test points(m)[end-1] == Point(0.95, 1, 1)
@test points(m)[end]   == Point(1.00, 1, 1)
