using Base.Test
using FAT.Utils
using HDF5


# use this as test case
casedir = "ldc_test"

# first overwrite the file
OpenFoamToHDF5(casedir; 
               fielddtype=Float64, 
               overwrite=true, 
               dimensions=2, 
               verbose=false)

# try to write it again
@test_throws ErrorException OpenFoamToHDF5(casedir)

# try different dimensions
@test_throws ErrorException OpenFoamToHDF5(casedir; 
                                           overwrite=true, 
                                           dimensions=1, 
                                           verbose=false)
@test_throws ErrorException OpenFoamToHDF5(casedir; 
                                           overwrite=true, 
                                           dimensions=4, 
                                           verbose=false)

# test type of saved data
for dtype in [Float32, Float64]
    OpenFoamToHDF5(casedir; 
                   fielddtype=dtype, 
                   overwrite=true, 
                   dimensions=2, 
                   verbose=false)
    fh = h5open("$casedir/data.h5")
    @test eltype(read(fh["1.000000/U/internalField"])) == dtype
    @test eltype(read(fh["1.000000/U/boundaryField"])) == dtype
    close(fh)
end

# check attribute for dimension
for D in ["2", "3"]
    OpenFoamToHDF5(casedir; 
                   overwrite=true, 
                   dimensions=parse(Int, D), verbose=false)
    fh = h5open("$casedir/data.h5")
    @test read(attrs(fh)["dimensions"]) == D
    close(fh)
end

# get times
times_string = FAT.Utils.get_snapshots_times(casedir)
times_float = [parse(Float64, t) for t in times_string]

# check everyone is in
for i in 0:1:10
    @test i in times_float
end

# check nothing else is in there
@test length(times_string) == length(0:10)

# check all times in the hdf5 file
fh = h5open("$casedir/data.h5")

# this will be a list of strings
h5_times = sort!([parse(Float64, i) for i in names(fh)])

# we test for equality, but remember that in 
# times_float we also have zero which we did not load
@test h5_times == times_float[2:end]

# close it
close(fh)

# now test the skip argument
OpenFoamToHDF5(casedir; 
               overwrite=true, 
               dimensions=2, 
               verbose=false, skip=2)

fh = h5open("$casedir/data.h5")
h5_times = sort!([parse(Float64, i) for i in names(fh)])
@test h5_times == times_float[3:end]

# close it
close(fh)

# skip cannot be too big
@test_throws  ErrorException OpenFoamToHDF5(casedir; 
               overwrite=true, 
               dimensions=2, 
               verbose=false, skip=-1)

@test_throws  ErrorException OpenFoamToHDF5(casedir; 
               overwrite=true, 
               dimensions=2, 
               verbose=false, skip=11)