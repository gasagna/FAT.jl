# ------------------------------------------------------------------- #
# Copyright 2015-2022, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

import LinearAlgebra
import Statistics

export AbstractField, 
       ScalarField, 
       VectorField, 
       TensorField, 
       grad, 
       grad!,
       curl, 
       curl!,
       mesh, 
       dotgrad, 
       dotgrad!, 
       projections


""" Abstract type for fields

    Type parameters
    ---------------
    D : number of spatial dimensions, either 2 or 3. This is useful
        because the gradient of a scalar field in D dimensions is a D-vector.
    T : type of the field values, e.g. Float64
    M : the type of Mesh

"""
abstract type AbstractField{D, T<:Real, M<:Mesh} <: AbstractVector{T} end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, u::AbstractField)
    print(io, "$(typeof(u)) object\n")
    print(io, "  ~ Spatial dimensions : $(ndims(u))\n")
    print(io, "  ~ Mesh information:\n")
    print(io, "    ~ ")
    show(io, mesh(u); gap="     ")
end


""" Type representing scalar fields, e.g. pressure, and velocity components.

    Notes
    -----
    Two dimensional simulations in OpenFOAM are inherently 3D, but on 
    special meshes. We could always have 3D data structure, but we would
    loose some efficiency in post-processing 2D simulations data.

    Fields
    ------
    internalField : a vector of the field values at the cell centres
    boundaryField : a vector, containing the values of boundary faces
             mesh : the mesh
"""
struct ScalarField{D, T, M} <: AbstractField{D, T, M}
    internalField::Vector{T}
    boundaryField::Vector{T}
    mesh::M
    function ScalarField(mesh::M, 
                            D::Integer, 
                internalField::AbstractVector{T}, 
                boundaryField::AbstractVector{T}) where {T, M}
        D in [2, 3] || error("`D` must be either 2 or 3")
        return new{D, T, M}(internalField, boundaryField, mesh)
    end
end


""" Type representing vector fields, e.g. velocity or vorticity.

    Notes
    -----
    ~ Type parameters D, T, M are as for the scalar field type.

    Fields
    ------
    scalars : a tuple of D ScalarField objects
       mesh : the mesh
"""
struct VectorField{D, T, M, S<:ScalarField{D, T, M}} <: AbstractField{D, T, M}
    scalars::NTuple{D, S}
    mesh::M
    function VectorField(mesh::M, scalars::NTuple{D, S}) where {D, T, M, S<:ScalarField{D, T, M}}
        return new{D, T, M, S}(scalars, mesh)
    end
end

# add syntax such as u[1] or u[:x]
@inline function Base.getindex(u::VectorField{D}, i::Union{Integer, Symbol}) where {D}
    D == 3 && (i == 3 || i == :w || i == :W || i == :z) && return u.scalars[3]
              (i == 2 || i == :v || i == :V || i == :y) && return u.scalars[2]
              (i == 1 || i == :u || i == :U || i == :x) && return u.scalars[1]
    throw(BoundsError("component `$i` not understood"))
end

""" Type representing tensor fields, e.g. the velocity gradient tensor

    Notes
    -----
    ~ Type parameters D, T, M are as for the scalar and vector fields types.

    Fields
    ------
    vectors : a tuple of D VectorField objects
       mesh : the mesh
"""
struct TensorField{D, T, M} <: AbstractField{D, T, M}
    vectors::NTuple{D, VectorField{D, T, M}}
    mesh::M
    function TensorField(mesh::Mesh, D::Integer, ::Type{T}=Float64) where {T<:Real}
        return new{D, T, typeof(mesh)}(ntuple(i->VectorField(mesh, D, T), D), mesh)
    end
end

# add syntax such as u[1] or u[:x]
@inline function Base.getindex(u::TensorField{D}, i::Union{Integer, Symbol}) where {D}
    D == 3 && (i == 3 || i == :w || i == :W || i == :z) && return u.vectors[3]
              (i == 2 || i == :v || i == :V || i == :y) && return u.vectors[2]
              (i == 1 || i == :u || i == :U || i == :x) && return u.vectors[1]
    throw(BoundsError("component `$i` not understood"))
end

# add syntax such as u[1] or u[:x]
@inline Base.getindex(u::TensorField, 
                      i::Union{Integer, Symbol},
                      j::Union{Integer, Symbol}) = u[i][j]

" Get the type of the field data "
Base.eltype(u::AbstractField{D, T}) where {D, T} = T

" Number of spatial dimensions "
Base.ndims(u::AbstractField{D}) where {D} = D

" Get the mesh of a field "
mesh(u::AbstractField) = u.mesh

# equality between fields
Base.:(==)(u::ScalarField, v::ScalarField) = 
    (u.internalField == v.internalField && u.boundaryField == v.boundaryField)
Base.:(==)(u::VectorField, v::VectorField) = u.scalars == v.scalars
Base.:(==)(u::TensorField, v::TensorField) = u.vectors == v.vectors

# similar
Base.similar(f::ScalarField) = ScalarField(f.mesh, ndims(f), similar(f.internalField), similar(f.boundaryField)) 
Base.similar(f::VectorField) = VectorField(f.mesh, similar.(f.scalars))

# zero and copy 
Base.zero(u::AbstractField) = (v = similar(u); v .= 0; v)
Base.copy(u::AbstractField) = (v = similar(u); v .= u; v)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Overload operators on ScalarFields ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~ BROADCASTING ~
const ScalarFieldStyle = Broadcast.ArrayStyle{ScalarField}
Base.BroadcastStyle(::Type{<:ScalarField}) = Broadcast.ArrayStyle{ScalarField}()

const VectorFieldStyle = Broadcast.ArrayStyle{VectorField}
Base.BroadcastStyle(::Type{<:VectorField}) = Broadcast.ArrayStyle{VectorField}()

# scalarfield
@inline function Base.materialize!(dest::ScalarField,
                                     bc::Broadcast.Broadcasted{ScalarFieldStyle})
    Base.materialize!(dest.internalField, unpack(bc, Val(:internalField)))
    Base.materialize!(dest.boundaryField, unpack(bc, Val(:boundaryField)))
    return dest
end

@inline unpack(x::ScalarField, ::Val{:internalField}) = x.internalField
@inline unpack(x::ScalarField, ::Val{:boundaryField}) = x.boundaryField

# vectorfield
@inline function Base.materialize!(dest::VectorField{D},
                                     bc::Broadcast.Broadcasted{VectorFieldStyle}) where {D}
    materialize!(dest[1], unpack(bc, Val(1)))
    materialize!(dest[2], unpack(bc, Val(2)))
    D == 3 && materialize!(dest[3], unpack(bc, Val(3)))
    return dest
end

@inline unpack(x::VectorField, ::Val{i}) where {i} = 
    (@inbounds v = x.scalars[i]; return v)

# generic
@inline unpack(bc::Broadcast.Broadcasted, 
             item::Val) = Broadcast.Broadcasted(bc.f, _unpack(bc.args, item))
@inline unpack(x::Any, ::Val) = x

@inline _unpack(args::Tuple, item) = (unpack(args[1], item), _unpack(Base.tail(args), item)...)
@inline _unpack(args::Tuple{Any}, item) = (unpack(args[1], item),)
@inline _unpack(args::Tuple{}, item) = ()