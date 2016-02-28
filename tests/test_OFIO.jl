using Base.Test
using FAT.OFIO

casedir = "/Users/davide/Dropbox/Codes/FAT/tests/ldc_test"

# MESH READERS
points_data = FAT.OFIO.reader(casedir, "points")
@test length(points_data) == 882
@test points_data[1] == (0.0, 0.0, 0.0)
@test points_data[end] == (1.0, 1.0, 1.0)

faces_data = FAT.OFIO.reader(casedir, "faces")
@test length(faces_data) == 1640
@test faces_data[1] == (1+1, 22+1, 463+1, 442+1)
@test faces_data[end] == (859+1, 860+1, 881+1, 880+1)

owners_data = FAT.OFIO.reader(casedir, "owner")
@test length(owners_data) == 1640
@test owners_data[1] == 0+1
@test owners_data[end] == 399+1

neighbour_data = FAT.OFIO.reader(casedir, "neighbour")
@test length(neighbour_data) == 760
@test neighbour_data[1] == 1+1
@test neighbour_data[end] == 399+1


# Some functions
@test FAT.OFIO.is_patch_name("hello") == true
@test FAT.OFIO.is_patch_name("   hello") == true
@test FAT.OFIO.is_patch_name("   hello you") == false
@test FAT.OFIO.is_patch_name("hello you") == false
@test FAT.OFIO.is_patch_name("hello you  ") == false

# read boundary file
out = FAT.OFIO.read_boundary(casedir)
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
# 	@test out[name][1] == FAT.Constants.BC_FIXEDVALUE
# 	@test out[name][2] == zeros(20, 2)
# end