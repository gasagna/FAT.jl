using Base.Test
using FAT.OFIO
using HVectors

let 
    casedir_b = "./ldc_test_binary"
    casedir_a = "./ldc_test"

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
end


# ~~~ helper functions ~~~
@test FAT.OFIO.fileformat("./ldc_test_binary/constant/" *
                          "polyMesh/points") == "binary"
@test FAT.OFIO.fileformat("./ldc_test/constant/" *
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

# ~~~~ Read internal vector field ~~~~

# @test FAT.OFIO.is_patch_name("hello") == true
# @test FAT.OFIO.is_patch_name("   hello") == true
# @test FAT.OFIO.is_patch_name("   hello you") == false
# @test FAT.OFIO.is_patch_name("hello you") == false
# @test FAT.OFIO.is_patch_name("hello you  ") == false

# # read boundary file
# out = FAT.OFIO.read_boundary(casedir)
# for t in [:top, :bottom, :left, :right, :front1, :back0]
#     @test t in keys(out)
# end
# # these are non-empty patches
# for t in [:top, :bottom, :left, :right]
#     @test out[t][1] == false
# end
# # these are empty patches
# for t in [:front1, :back0]
#     @test out[t][1] == true
# end


# # DATA FILE READERS
# # data reader for pressure
# fh = open(joinpath(casedir, "1/p"), "r")
# p = FAT.OFIO.read_internal_scalar_field(fh, Float64)
# @test p[1] == 2.27405e-5
# @test p[end] == 46.4739
# close(fh)

# # data reader for velocity
# fh = open(joinpath(casedir, "1/U"), "r")
# # read first two velocity components
# U = FAT.OFIO.read_internal_vector_field(fh, 2, Float64)
# @test vec(U[1, :]) == [0.000259015, -0.000255147]
# # fixes bug when exponential notation is used 
# @test vec(U[21, :]) == [-6.3677e-5, -0.000192864]
# @test vec(U[end, :]) == [0.30517, -0.147796]

# # read three velocity components
# seek(fh, 0)
# U = FAT.OFIO.read_internal_vector_field(fh, 3, Float64)
# @test vec(U[1, :]) == [0.000259015, -0.000255147, 0.0]
# close(fh)

# # read boundary field
# fh = open(joinpath(casedir, "1/U"), "r")
# out = FAT.OFIO.read_boundary_vector_field(casedir, fh, 2, Float64)
# println(out)

# @test out[:top][1] == FAT.Constants.BC_FIXEDVALUE
# @test out[:top][2][:, 1] == ones(20)
# @test out[:top][2][:, 2] == zeros(20)

# for name in [:left, :right, :bottom]
#   @test out[name][1] == FAT.Constants.BC_FIXEDVALUE
#   @test out[name][2] == zeros(20, 2)
# end