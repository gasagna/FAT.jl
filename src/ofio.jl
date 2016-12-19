# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module OFIO

using HeterogeneousVectors


# ~~~~~~~~~~~~~~~~~
# Helper functions
# ~~~~~~~~~~~~~~~~~

""" 
        iscasedir(casedir)

    Check if `directory` is an OpenFOAM case directory.
"""
function iscasedir(casedir::AbstractString)
    files = readdir(casedir)
    "constant" in files && "system" in files
end

""" 
        gotomatch(f, regex)

    Scan through an open file `f` until a line match 
    with `regex`. Then do nothing. It is an error if 
    no line matches.
"""
gotomatch(f::IO, regex::Regex) = (matchingline(f, regex); nothing)

""" 
        matchingline(f, regex)

    Scan through an open file `f` until a line match 
    with `regex`. Then return the matching line. It is 
    an error if no line matches.
"""
function matchingline(f::IO, regex::Regex) 
    for line in eachline(f)
        ismatch(regex, line) && return line
    end
    error("no match found with $regex in $f")
end

""" 
        fileformat(filename)

    Detect file format of an OpenFOAM file by parsing the 
    FoamFile dictionary at the beginning of the the file.
    An error is raised if no format is detected.
"""
function fileformat(filename::AbstractString)
    open(filename) do f
        m = matchall(r"binary|ascii", matchingline(f, r"format"))
        length(m) == 0 && error("no suitable format found in $filename")
        return m[1]
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to read OpenFOAM mesh files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The below functions are specific implementations
# for each type of file ("points", "faces"), and for 
# binary/ascii format.
   
# ~~~ binary OpenFOAM mesh readers

function read_points_binary(casedir::AbstractString, mtype::Type)
    open(joinpath(casedir, "constant/polyMesh/points")) do f
        # ~~~ allocate times three as we have 3d points ~~~
        N = parse(Int, matchingline(f, r"^[0-9]+"))
        out = Vector{Float64}(3*N)
        read(f, Char); read!(f, out)
        
        # ~~~ we need a vector of 3-tuples ~~~
        reinterpret(NTuple{3, mtype}, out)
    end
end

function read_faces_binary(casedir::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/faces")) do f
        # ~~~ read indices ~~~
        N = parse(Int, matchingline(f, r"^[0-9]+"))
        idxs = Vector{UInt32}(N)
        read(f, Char); read!(f, idxs)

        # ~~~ read data ~~~
        # There are two lists in this file
        N = parse(Int, matchingline(f, r"^[0-9]+"))
        data = Vector{UInt32}(N)
        read(f, Char); read!(f, data)

        # ~~~ add one to data ~~~
        @simd for i in eachindex(data)
            @inbounds data[i] += 1
        end
        # ~~~ return heterogeneous vector ~~~
        HVector(data, idxs)
    end
end

# generic function valid for "owner" and "neighbour"
function read_on_binary(casedir::AbstractString, fname::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/", fname)) do f
        # ~~~ read data  ~~~
        N = parse(Int, matchingline(f, r"^[0-9]+"))
        out = Vector{UInt32}(N)
        read(f, Char); read!(f, out)

        # ~~~ add one to data ~~~
        @simd for i in eachindex(out) 
            @inbounds out[i] += 1
        end
        out
    end
end


# ~~~ ASCII OpenFOAM mesh readers

function read_points_ascii(casedir::AbstractString, mtype::Type)
    open(joinpath(casedir, "constant/polyMesh/points")) do f
        # ~~~ number of lines to read ~~~
        N = parse(Int, matchingline(f, r"^[0-9]+"))
        # skip (
        readline(f); 
        # preallocate
        out = Vector{NTuple{3, mtype}}(N)

        # ~~~ parse all rows ~~~
        for i = 1:length(out)
            m = matchall(r"-?[\d.]+(?:e-?\d+)?", readline(f))
            @inbounds out[i] = (parse(mtype, m[1]), 
                                parse(mtype, m[2]), 
                                parse(mtype, m[3]))
        end
        out
    end
end

function read_faces_ascii(casedir::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/faces")) do f
        # ~~~ number of lines to read ~~~
        N = parse(Int, matchingline(f, r"^[0-9]+"))
        # skip (
        readline(f); 

        # faces can have any number of points, so we use an HVector
        out = HVector{UInt32}()

        # ~~~ parse all rows ~~~
        for i = 1:N
            # parse line: this is the most intensive bit
            m = split(readline(f), [' ', '(', ')'],  keep=false)
            # Loop over number of face points. We add 0x00000001 
            # because we want a 1-based data structure
            @unlocked out for i = 1:parse(Int, m[1])
                push!(out, parse(UInt32, m[i+1]) + 0x00000001)
            end
        end
        out
    end
end

# generic function valid for "owner" and "neighbour"
function read_on_ascii(casedir::AbstractString, filename::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/", filename)) do f
        # ~~~ number of lines to read ~~~
        N = parse(Int, matchingline(f, r"^[0-9]+"))
        # skip (
        readline(f); 
        # ~~~ parse all rows ~~~
        UInt32[parse(UInt32, readline(f)) + one(UInt32) for i = 1:N]
    end
end

# also define methods such as read_owner_binary, ...
for fmt in [:_ascii, :_binary], filename in ["owner", "neighbour"]
    @eval $(symbol(:read_, filename, fmt))(casedir) = 
        $(symbol(:read_on, fmt))(casedir, $filename)
end

# Generic functions that are called in mesh.jl and does dispatch, 
# depending on the type of file and file format.

"""
        reader(casedir, filename, args...)

    Generic entry-point function for reading OpenFOAM mesh files. 
    This function is used in meshes.jl for loading all types of 
    mesh files. The last argument args... can be empty or have some
    value, as the "_points" functions expect an "mtype" argument, 
    whereas the other do not.
"""
function reader(casedir::AbstractString, filename::AbstractString, args...)
    fmt = fileformat(joinpath(casedir, "constant/polyMesh/$filename"))
    fname = symbol("read_", filename, "_", fmt) # define function name
    @eval $fname($casedir, $args...)
end
      
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to read patches
# ~~~~~~~~~~~~~~~~~~~~~~~~~
""" A line contains a patch name only if it contains a single word. 

    See http://stackoverflow.com/questions/8536749/regex-to
        -check-a-string-contains-just-one-word
"""
is_patch_name(line::AbstractString) = length(matchall(r"\w+", line)) == 1

""" Parse the `boundary` file and return a dictionary structure.

    Notes
    -----
    The dictionary is indexed by symbols, the patch names. The values are 
    (Bool, UInt32, UInt32) tuples. The boolean value tells whether the face
    is an empty face. The two UInts are the number of faces in the patch 
    and the id of the starting face, with 1-based indexing.
"""
function read_boundary(casedir::AbstractString)
    # initialise dict
    patches = Dict{Symbol, Tuple{Bool, UInt32, UInt32}}()
    # open file
    f = open(joinpath(casedir, "constant/polyMesh/boundary"), "r")
    # go to line where number of patches is shown
    npatches = parse(Int, matchingline(f, r"^[0-9]+"))
    while !eof(f)
        line = readline(f)
        if is_patch_name(line)
            patchname = symbol(strip(line))
            isempty = contains(matchingline(f, r"type"), "empty")
            nfaces = parse(UInt32, split(strip(matchingline(f, r"nFaces"), [' ', ';', '\n']))[2])
            startface = parse(UInt32, split(strip(matchingline(f, r"startFace"), [' ', ';', '\n']))[2])
            patches[patchname] = (isempty, nfaces, startface+UInt32(1))
        end
    end 
    patches
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to read OpenFOAM vector field files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""" Read an OpenFoam vector field file and return the internal 
    and boundary fields, as vectors of vectors of appropriate size.

    `dims` is a tuple with the dimensions that will be read. 
"""
function read_vector_field(filename::AbstractString, dims::Tuple{Vararg{Int}})
    open(filename, "r") do f
        # Parse file format. Raises an error if not ascii or binary
        format = fileformat(f) 

        # Read internal and boundary fields. We assume that the former
        # preceeds the latter in the output files. Otherwise we fail 
        # badly. We could do dispatch based on the format, but we
        # would need to do dynamic dispatch at runtime. The if clause
        # here leads to faster code, although it does not really matter
        # for the overall performance.
        format == "binary" && 
            return (read_internal_vector_field_binary(f, dims),
                    read_boundary_vector_field_binary(f, dims))
        
        format == "ascii" && 
            return (read_internal_vector_field_ascii(f, dims),
                    read_boundary_vector_field_ascii(f, dims))

    end
end

function read_internal_vector_field_ascii(f::IO, dims::Tuple{Vararg{Int}})
    gotomatch(f, r"internalField")
    nlines = parse(Int, readline(f)); readline(f)
    # load data into a large matrix
    out = zeros(nlines, length(dims))
    for i = 1:nlines
        m = matchall(r"-?[\d.]+(?:e-?\d+)?", readline(f))
        for j = 1:length(dims)
            out[i, j] = parse(Float64, m[dims[j]])
        end
    end
    [slice(out, :, j) for j in 1:length(dims)]
end

function read_internal_vector_field_binary(f::IO, dims::Tuple{Vararg{Int}})
    gotomatch(f, r"internalField")
    # get number of internal cell centers
    ncells = parse(Int, readline(f))
    # skip the leading (
    read(f, Char)
    # read into large matrix, then split
    # TODO: have this more efficient
    out = read(f, Float64, (3, ncells))
    [slice(out, d, :) for d in dims]
end

function read_boundary_vector_field_ascii(casedir::AbstractString, f::IO, dims::Tuple{Vararg{Int}})
    # Read the boundary file. This happens every time!
    patches = read_boundary(casedir)
    # value[2] is the number of faces, so we sum them
    nboundaryfaces = sum([Int(value[2]) for (key, value) in patches])
    # number of internal faces. value[3] is the starting face, so we get the min
    nInternalFaces = minimum([Int(value[3]) for (key, value) in patches]) - 1
    output = zeros(nboundaryfaces, length(dims))
    gotomatch(f, r"boundaryField")
    while !eof(f)
        line = readline(f)
        # go to a patch name
        if is_patch_name(line)
            patchname = symbol(strip(line))
            _, patchtype = split(strip(matchingline(f, r"type"), [' ', ';', '\n']))
            # number of faces in the patch. See the `read_boundary` function above.
            n = Int(patches[patchname][2])
            # for empty patches we do nothing
            if patchtype == "empty"
                nothing 
            # the easiest case to deal with, when velocity has a fixed value
            elseif patchtype == "fixedValue"
                # go to next line and parse the () part. Parse all components
                g = match(r"\(([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+) ([-+]?[0-9]*\.?[0-9]+)\)", readline(f))
                # set all entries to the same value 
                val =  [parse(Float64, g[d]) for d in dims]
                for i = patches[patchname][3]:(patches[patchname][3]+patches[patchname][2]-1)
                    # we need to remove the number of internal faces from i
                    i_ = i - nInternalFaces
                    output[i_, :] = val
                end
            elseif patchtype == "calculated"
                # read number of lines that need to be read. Note that the ( is skipped
                nlines = goToGoodLine(f)
                for i = patches[patchname][3]:(patches[patchname][3]+patches[patchname][2]-1)
                    m = matchall(r"-?[\d.]+(?:e-?\d+)?", readline(f))
                    # we need to remove the number of internal faces from i
                    i_ = i - nInternalFaces
                    for j = 1:length(dims)
                        output[i_, j] = parse(Float64, m[dims[j]])
                    end
                end
            else
                # if we do not know the type of the patch we raise an error
                #error("patchtype $patchtype not understood")
                nothing
            end
        end
    end 
    [slice(output, :, j) for j in 1:length(dims)]
end

function read_boundary_vector_field_binary(casedir::AbstractString, f::IO, dimensions::Integer)
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to read OpenFOAM scalar field files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function read_scalar_field(filename::AbstractString)
end

function read_internal_scalar_field_ascii(f::IO)
    gotomatch(f, r"internalField")
    nlines = parse(Int, readline(f)); readline(f)
    mtype[parse(mtype, readline(f)) for i in 1:nlines]
end

function read_internal_scalar_field_binary(f::IO)
end


end