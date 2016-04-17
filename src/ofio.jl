# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module OFIO

using HVectors


# ~~~~~~~~~~~~~~~~~
# Generic functions
# ~~~~~~~~~~~~~~~~~

""" 
        iscasedir(casedir)

    Check if `directory` is an OpenFoam case directory.
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
gotomatch(f::IOStream, regex::Regex) = (linematches(f, regex); nothing)

""" 
        linematches(f, regex)

    Scan through an open file `f` until a line match 
    with `regex`. Then return the matching line. It is 
    an error if no line matches.
"""
function linematches(f::IOStream, regex::Regex) 
    for line in eachline(f)
        ismatch(regex, line) && return line
    end
    error("no match found with $regex in $f")
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to read OpenFoam mesh files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
# binary OpenFOAM mesh reader

""" Read the `points` file in the constant/polyMesh directory """
function read_points_binary(casedir::AbstractString, mtype::Type)
    open(joinpath(casedir, "constant/polyMesh/points")) do f
        # ~~~ allocate times three as we have 3d points ~~~
        N = parse(Int, linematches(f, r"^[0-9]+"))
        out = Vector{Float64}(3*N)
        read(f, Char); read!(f, out)
        
        # ~~~ we need a vector of 3-tuples ~~~
        reinterpret(NTuple{3, mtype}, out)
    end
end

""" Read the `faces` file in the constant/polyMesh directory """
function read_faces_binary(casedir::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/faces")) do f
        # ~~~ read indices ~~~
        N = parse(Int, linematches(f, r"^[0-9]+"))
        idxs = Vector{UInt32}(N)
        read(f, Char); read!(f, idxs)

        # ~~~ read data ~~~
        # There are two lists in this file
        N = parse(Int, linematches(f, r"^[0-9]+"))
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

""" Read the `owner` or `neighbour` file in the constant/polyMesh directory """
function read_on_binary(casedir::AbstractString, fname::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/", fname)) do f
        # ~~~ read data  ~~~
        N = parse(Int, linematches(f, r"^[0-9]+"))
        out = Vector{UInt32}(N)
        read(f, Char); read!(f, out)

        # ~~~ add one to data ~~~
        @simd for i in eachindex(out) 
            @inbounds out[i] += 1
        end
        out
    end
end


# ASCII OpenFOAM mesh reader

""" Read the `points` file in the constant/polyMesh directory """
function read_points_ascii(casedir::AbstractString, mtype::Type)
    open(joinpath(casedir, "constant/polyMesh/points")) do f
        # ~~~ number of lines to read ~~~
        N = parse(Int, linematches(f, r"^[0-9]+"))
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

""" Read the `faces` file in the constant/polyMesh directory """
function read_faces_ascii(casedir::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/faces")) do f
        # ~~~ number of lines to read ~~~
        N = parse(Int, linematches(f, r"^[0-9]+"))
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

""" Read the `owner` or `neighbour` file in the constant/polyMesh directory """
function read_on_ascii(casedir::AbstractString, filename::AbstractString)
    open(joinpath(casedir, "constant/polyMesh/", filename)) do f
        # ~~~ number of lines to read ~~~
        N = parse(Int, linematches(f, r"^[0-9]+"))
        # skip (
        readline(f); 
        # ~~~ parse all rows ~~~
        UInt32[parse(UInt32, readline(f)) + one(UInt32) for i = 1:N]
    end
end


# Generic functions that are called in mesh.jl and does dispatch, 
# depending on the type of file, i.e. ascii or binary

function read_points(casedir::AbstractString, mtype::Type)
    f = open(joinpath(casedir, "constant/polyMesh/points"))
    line = linematches(f, r"format") 
    close(f)
    contains(line, "binary") && return read_points_binary(casedir, mtype)
    contains(line, "ascii")  && return read_points_ascii(casedir, mtype)
    error("no suitable format found [ascii/binary]")
end

function read_faces(casedir::AbstractString)
    f = open(joinpath(casedir, "constant/polyMesh/faces"))
    line = linematches(f, r"format") 
    close(f)
    contains(line, "binary") && return read_faces_binary(casedir)
    contains(line, "ascii")  && return read_faces_ascii(casedir)
    error("no suitable format found [ascii/binary]")
end

function read_on(casedir::AbstractString, filename::AbstractString)
    f = open(joinpath(casedir, "constant/polyMesh/$filename"))
    line = linematches(f, r"format") 
    close(f)
    contains(line, "binary") && return read_on_binary(casedir, filename)
    contains(line, "ascii")  && return read_on_ascii(casedir, filename)
    error("no suitable format found [ascii/binary]")
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
    npatches = parse(Int, linematches(f, r"^[0-9]+"))
    while !eof(f)
        line = readline(f)
        if is_patch_name(line)
            patchname = symbol(strip(line))
            isempty = contains(linematches(f, r"type"), "empty")
            nfaces = parse(UInt32, split(strip(linematches(f, r"nFaces"), [' ', ';', '\n']))[2])
            startface = parse(UInt32, split(strip(linematches(f, r"startFace"), [' ', ';', '\n']))[2])
            patches[patchname] = (isempty, nfaces, startface+UInt32(1))
        end
    end 
    patches
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to read OpenFoam output files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


""" Read scalar internal field """
function read_internal_scalar_field(f::IOStream, mtype::Type=Float64)
    gotomatch(f, r"internalField")
    nlines = parse(Int, readline(f)); readline(f)
    mtype[parse(mtype, readline(f)) for i in 1:nlines]
end

""" Read vector internal field.

    Notes
    -----
    The argument `dimensions` specifies how many components will be read.

"""
function read_internal_vector_field(f::IOStream, dimensions::Integer, mtype::Type=Float64)
    gotomatch(f, r"internalField")
    nlines = parse(Int, readline(f)); readline(f)
    out = Matrix{mtype}(nlines, dimensions)
    for i = 1:nlines
        m = matchall(r"-?[\d.]+(?:e-?\d+)?", readline(f))
        for j = 1:dimensions
            out[i, j] = parse(mtype, m[j])
        end
    end
    out
end

""" Read the vector field on the domain boundary.

    Currently we recognize the following types of values:

    ~ type fixedValue
      value uniform (ux, uy, uz)

    We only support clearly formatted files, so do not mess 
    up with the output files.
"""
function read_boundary_vector_field(casedir::AbstractString, f::IOStream, dimensions::Integer, mtype::Type=Float64)
    # Read the boundary file. This happens every time!
    patches = read_boundary(casedir)
    # value[2] is the number of faces, so we sum them
    nboundaryfaces = sum([Int(value[2]) for (key, value) in patches])
    # number of internal faces. value[3] is the starting face, so we get the min
    nInternalFaces = minimum([Int(value[3]) for (key, value) in patches]) - 1
    output = zeros(mtype, nboundaryfaces, dimensions)
    gotomatch(f, r"boundaryField")
    while !eof(f)
        line = readline(f)
        # go to a patch name
        if is_patch_name(line)
            patchname = symbol(strip(line))
            _, patchtype = split(strip(linematches(f, r"type"), [' ', ';', '\n']))
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
                val =  [parse(mtype, g[i]) for i in 1:dimensions]
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
                    for j = 1:dimensions
                        output[i_, j] = parse(mtype, m[j])
                    end
                end
            else
                # if we do not know the type of the patch we raise an error
                #error("patchtype $patchtype not understood")
                nothing
            end
        end
    end 
    output
end

end