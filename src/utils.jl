# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Utils

export OpenFoamToHDF5

import HDF5: h5open, 
			 g_create, 
			 close, 
			 attrs

import FAT.OFIO: read_internal_vector_field, 
				 read_boundary_vector_field

""" Create an array of time strings at which a snapshot is available """
get_snapshots_times(casedir::AbstractString) = 
	sort(filter(f->ismatch(r"^[0-9]*\.?[0-9]+$", f), 
		 readdir(casedir)), by=float) 

""" Read OpenFoam output and convert it to a single large HDF5 
	file in the `casedir` directory.

	Only the velocity vector field is loaded, as pressure should not
	be so useful, unless one wants to calculate forces.

	Parameters
	----------
	   casedir : the OpenFoam case directory
	fielddtype : type of the floating point data used for storage
	 overwrite : whether to overwrite the data file if it exists
	dimensions : either 2 or 3. The number of useful spatial components
				 of the simulation. If 3, all three velocity components
				 are read and saved. If 2 only the first two, i.e. those 
				 along the x and y directions. Defaults to 2.
	   verbose : if true, print current time being read
"""
function OpenFoamToHDF5(casedir::AbstractString; 
						toload::Vector{Symbol}=[:U, :W],
						skip::Integer=1,
					    fielddtype::Type=Float64, 
					    overwrite::Bool=false, 
					    dimensions::Integer=2, 
					    verbose=true)
	# FIXME: check if it is better to have multiple hdf files
	# check dimensions
	dimensions in [2, 3] || error("`dimensions` must be either 2 or 3")

	# make filename
	filename = joinpath(casedir, "data.h5")

	# check not to overwrite
	overwrite || ispath(filename) && error("$filename already exists.")
	
	# open HDF5 main file and set number of dimensions in an attribute
	fh = h5open(filename, "w")
	attrs(fh)["dimensions"] = string(dimensions)
	try 
		# get times
		times = get_snapshots_times(casedir)
		0 <= skip < length(times) || error("value of `skip` not meaningful")

		for (i, time) in enumerate(times)
			# skip the first `skip` directories, as we might want to skip the 
			# initial condition, whose files are difficult to parse and we 
			# might not care
			i <= skip && continue

			# use six decimal digits for the groups names
			ts = @sprintf "%.6f" parse(Float64, time)

			# debug output
			verbose && print("\r Reading t=$(ts)"); flush(STDOUT)
			
			# create groups for each time
			time_group = g_create(fh, ts)

			# we need to remember how big is the boundary file
			N = 0 
			if :U in toload 
				U_group = g_create(time_group, "U")
				U_group["internalField"] = read_internal_vector_field(
						open(joinpath(casedir, time, "U")), dimensions, fielddtype)
				data = read_boundary_vector_field(
						casedir, open(joinpath(casedir, time, "U")), 
						dimensions, fielddtype)
				N = length(data)
				U_group["boundaryField"] = data
			end
			if :W in toload 
				W_group = g_create(time_group, "W")
				W_group["internalField"] = read_internal_vector_field(
						open(joinpath(casedir, time, "W")), dimensions, fielddtype)
				W_group["boundaryField"] = 0.0*zeros(fielddtype, N, dimensions)
			end
			if :vorticity in toload 
				U_group = g_create(time_group, "vorticity")
				# For 2D flows we only need the last component
				a = read_internal_vector_field(
					open(joinpath(casedir, time, "vorticity")), 3, fielddtype)
				U_group["internalField"] = dimensions == 2 ? a[:, end] : a
				b = read_boundary_vector_field(
					casedir, open(joinpath(casedir, time, "vorticity")), 3, 
					fielddtype)
				U_group["boundaryField"] = dimensions == 2 ? b[:, end] : b
			end
		end
		verbose && println("\n")
	finally
		close(fh)
	end
	nothing
end

end


