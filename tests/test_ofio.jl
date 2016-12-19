using Base.Test
using FAT.OFIO
using FAT.Meshes
using HeterogeneousVectors

casedir_b = "./ldc_test_binary"
casedir_a = "./ldc_test_ascii"

pts_a = FAT.OFIO.reader(casedir_a, "points", Float64)
pts_b = FAT.OFIO.reader(casedir_b, "points", Float64)
@test pts_a == pts_b

@test typeof(pts_a) == Vector{NTuple{3, Float64}}
@test length(pts_a) == 882
@test pts_a[1] == (0.0, 0.0, 0.0)
@test pts_a[end] == (1.0, 1.0, 1.0)

faces_a = FAT.OFIO.reader(casedir_a, "faces")
faces_b = FAT.OFIO.reader(casedir_b, "faces")
@test faces_a == faces_b

@test typeof(faces_a) == HVector{UInt32, UInt32}
@test length(faces_a) == 1640
@test faces_a[1] == [1+1, 22+1, 463+1, 442+1]
@test faces_a[end] == [859+1, 860+1, 881+1, 880+1]

owner_a = FAT.OFIO.reader(casedir_a, "owner")
owner_b = FAT.OFIO.reader(casedir_b, "owner")
@test owner_a == owner_b

@test typeof(owner_a) == Vector{UInt32}
@test length(owner_a) == 1640
@test owner_a[1] == 0+1
@test owner_a[end] == 399+1

neighbour_a = FAT.OFIO.reader(casedir_a, "neighbour")
neighbour_b = FAT.OFIO.reader(casedir_b, "neighbour")
@test neighbour_a == neighbour_b

@test typeof(neighbour_a) == Vector{UInt32}
@test length(neighbour_a) == 760
@test neighbour_a[1] == 1+1
@test neighbour_a[end] == 399+1


# ~~~ helper functions ~~~
@test FAT.OFIO.fileformat("./ldc_test_binary/constant/" *
                          "polyMesh/points") == "binary"
@test FAT.OFIO.fileformat("./ldc_test_ascii/constant/" *
                          "polyMesh/points") == "ascii"
@test_throws ErrorException FAT.OFIO.fileformat("./noformat")

@test FAT.OFIO.iscasedir("./ldc_test_binary") == true
@test FAT.OFIO.iscasedir("../") == false

# matchingline 
io = IOBuffer(
"""
hello
2
(
(0 0 0)
(0.05 0 0)
)

3
(
(0 0 0)
(0 0 0)
(0 0 0)
)
""")

@test FAT.OFIO.matchingline(io, r"^[0-9]+") == "2\n"
seek(io, 0)

@test FAT.OFIO.matchingline(io, r"hello") == "hello\n"
seek(io, 0)

@test_throws ErrorException FAT.OFIO.matchingline(io, r"not-matching")
seek(io, 0)

@test FAT.OFIO.matchingline(io, r"^[0-9]+") == "2\n"
@test FAT.OFIO.matchingline(io, r"^[0-9]+") == "3\n"
seek(io, 0)

# gotomatch
@test FAT.OFIO.gotomatch(io, r"hello") == nothing
@test readline(io) == "2\n"
seek(io, 0)

@test_throws ErrorException FAT.OFIO.gotomatch(io, r"not-matching")
seek(io, 0)


@test FAT.OFIO.is_patch_name("hello") == true
@test FAT.OFIO.is_patch_name("   hello") == true
@test FAT.OFIO.is_patch_name("   hello you") == false
@test FAT.OFIO.is_patch_name("hello you") == false
@test FAT.OFIO.is_patch_name("hello you  ") == false

# read boundary file
out = FAT.OFIO.read_boundary(casedir_a)
for t in [:top, :bottom, :left, :right, :front1, :back0]
    @test t in keys(out)
end
# these are non-empty patches
for t in [:top, :bottom, :left, :right]
    @test out[t][1] == false
end
# these are empty patches
for t in [:front1, :back0]
    @test out[t][1] == true
end

# ~~~~ Read internal vector field ~~~~

# compare between ascii and binary
for dims in [(1, 2), (2, 3), (1, 3), (1, 2, 3)]
    # open fields
    fh_a = open(joinpath(casedir_a, "1/U"), "r")
    fh_b = open(joinpath(casedir_b, "1/U"), "r")
    # load velocity components
    U_a = FAT.OFIO.read_internal_vector_field_ascii(fh_a, dims)
    U_b = FAT.OFIO.read_internal_vector_field_binary(fh_b, dims)

    for i = 1:length(dims)
        @test U_a[i] == U_b[i]
    end
end

# test individual values on ascii file
fh_a = open(joinpath(casedir_a, "1/U"), "r")
U_a = FAT.OFIO.read_internal_vector_field_ascii(fh_a, (1, 2, 3))
@test U_a[1][1:2] == [ 0.000259015397, 0.000163142844]
@test U_a[2][1:2] == [-0.000255146514, 0.000103467574]
@test U_a[3][1:2] == [ 0.000000000000, 0.000000000000]

# fixes bug when exponential notation is used 
@test U_a[1][21] == -6.36769612e-05


# ~~~~ Read boundary vector field ~~~~

# test length of output is consistent with input dimensions
for dims in [(1, 2), (2, 3), (1, 3), (1, 2, 3)]
    open(joinpath(casedir_a, "1/U"), "r") do fh_a
        out = FAT.OFIO.read_boundary_vector_field_ascii(casedir_a, fh_a, dims)
        @test length(out) == length(dims)
    end
end

# test individual values
msh = Mesh(casedir_a)
open(joinpath(casedir_a, "1/U"), "r") do fh_a
    out = FAT.OFIO.read_boundary_vector_field_ascii(casedir_a, fh_a, (1, 2, 3))
    # unitary velocity of upper wall 
    utop = ones(nfaces(patches(msh)[:top]))
    # zero velocity on all others
    uother = zeros(nboundaryfaces(msh) - nfaces(patches(msh)[:top]))
    @test out[1] == vcat(utop, uother)
    @test out[2] == zeros(nboundaryfaces(msh))
    @test out[3] == zeros(nboundaryfaces(msh))
end