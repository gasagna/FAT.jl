# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Fields

import Base: call, 
             ==,
             *, 
             /, 
             +,
             -,
             getindex,
             norm,
             show,
             eltype,
             ndims,
             zero,
             mean,
             fill!,
             similar

using FAT.Meshes

export AbstractField, 
       ScalarField, 
       VectorField, 
       TensorField, 
       inner, 
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
abstract AbstractField{D, T<:Real, M<:Mesh}

function show(io::IO, u::AbstractField)
    print(io, "$(typeof(u)) object at $(object_id(u))\n")
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
type ScalarField{D, T, M} <: AbstractField{D, T, M}
    internalField::Vector{T}
    boundaryField::Vector{T}
    mesh::M
    # the problem is that there is no way to infer the parameter D from 
    # the input arguments so D must be always specified. FIXME. WONTFIX
end

function ScalarField{T<:Real}(mesh::Mesh, D::Integer, dtype::Type{T}=Float64)
    D in [2, 3] || error("`D` must be either 2 or 3")
    ScalarField{D, dtype, typeof(mesh)}(Vector{dtype}(ncells(mesh)), 
                                        Vector{dtype}(nboundaryfaces(mesh)), 
                                        mesh)
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
type VectorField{D, T, M} <: AbstractField{D, T, M}
    scalars::NTuple{D, ScalarField{D, T, M}}
    mesh::M
end

VectorField{T<:Real}(mesh::Mesh, D::Integer, dtype::Type{T}=Float64) =
    VectorField(ntuple(i->ScalarField(mesh, D, dtype), D), mesh)

# add syntax such as u[1] or u[:x]
@inline function getindex{D}(u::VectorField{D}, i::Union{Integer, Symbol}) 
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
type TensorField{D, T, M} <: AbstractField{D, T, M}
    vectors::NTuple{D, VectorField{D, T, M}}
    mesh::M
end

TensorField{T<:Real}(mesh::Mesh, D::Integer, dtype::Type{T}=Float64) =
    TensorField(ntuple(i->VectorField(mesh, D, dtype), D), mesh)

# add syntax such as u[1] or u[:x]
@inline function getindex{D}(u::TensorField{D}, i::Union{Integer, Symbol}) 
    D == 3 && (i == 3 || i == :w || i == :W || i == :z) && return u.vectors[3]
              (i == 2 || i == :v || i == :V || i == :y) && return u.vectors[2]
              (i == 1 || i == :u || i == :U || i == :x) && return u.vectors[1]
    throw(BoundsError("component `$i` not understood"))
end

# add syntax such as u[1] or u[:x]
@inline getindex{D}(u::TensorField{D}, 
                    i::Union{Integer, Symbol},
                    j::Union{Integer, Symbol}) = u[i][j]

" Get the type of the field data "
eltype{D, T}(u::AbstractField{D, T}) = T

" Number of spatial dimensions "
ndims{D}(u::AbstractField{D}) = D

" Get the mesh of a field "
mesh(u::AbstractField) = u.mesh

# equality between fields
==(u::ScalarField, v::ScalarField) = (u.internalField == v.internalField && 
                                      u.boundaryField == v.boundaryField)
==(u::VectorField, v::VectorField) = u.scalars == v.scalars
==(u::TensorField, v::TensorField) = u.vectors == v.vectors

# similar, fill! and zero
similar{D, T}(u::ScalarField{D, T}) = ScalarField(mesh(u), D, T)
similar{D, T}(u::VectorField{D, T}) = VectorField(mesh(u), D, T)
similar{D, T}(u::TensorField{D, T}) = TensorField(mesh(u), D, T)

function fill!{D, T}(u::ScalarField{D, T}, val::Real) 
    fill!(u.internalField, val); fill!(u.boundaryField, val)
    u
end
zero{D, T}(u::ScalarField{D, T}) = fill!(ScalarField(mesh(u), D, T), zero(T))

function fill!{D, T}(u::VectorField{D, T}, val::Real)
    for d = 1:D
        fill!(u[d], zero(T))
    end
    u
end
zero{D, T}(u::VectorField{D, T}) = fill!(VectorField(mesh(u), D, T), zero(T))

function fill!{D, T}(u::TensorField{D, T}, val::Real)
    for d = 1:D
        fill!(u[d], zero(T))
    end
    u
end
zero{D, T}(u::TensorField{D, T}) = fill!(TensorField(mesh(u), D, T), zero(T))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Overload operators on ScalarFields ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Notes: 
#  ~ we never check if dimensions and mesh of the operands match
#  ~ only a subset of all possible operations have been overloaded
#    i.e., only those that were found useful at some point. Others 
#    can be added as well.

# scalar - scalar operations
for (op, fname) in zip([:-, :+, :*], [:sub!, :add!, :mul!])
    # fast in-place version scalar-scalar
    @eval function $fname(u::ScalarField, v::ScalarField, out::ScalarField)
              a = out.internalField; b = u.internalField; c = v.internalField  
              @simd for i in eachindex(a)
                        @inbounds a[i] = $op(b[i], c[i])
                    end
              a = out.boundaryField; b = u.boundaryField; c = v.boundaryField
              @simd for i in eachindex(a) 
                        @inbounds a[i] = $op(b[i], c[i])
                    end
              out
          end 
    # memory-allocating versions
    @eval $op(u::ScalarField, v::ScalarField) = $fname(u, v, similar(u))
end   

# Compute out = u*v + w. This is used in the 
# `dotgrad` function.
function muladd!(u::ScalarField,   v::ScalarField, 
                 w::ScalarField, out::ScalarField)
    a = out.internalField; b = u.internalField
    c = v.internalField;   d = w.internalField  
    @simd for i in eachindex(a)
        @inbounds a[i] = b[i]*c[i] + d[i]
    end
    a = out.boundaryField; b = u.boundaryField
    c = v.boundaryField;   d = w.boundaryField  
    @simd for i in eachindex(a)
        @inbounds a[i] = b[i]*c[i] + d[i]
    end
    out
end   

# scalar - real operations
for (op, fname) in zip([:*, :/], [:mul!, :div!])    
    @eval function $fname(u::ScalarField, v::Real, out::ScalarField)
          a = out.internalField; b = u.internalField
          @simd for i in eachindex(a)
                    @inbounds a[i] = $op(b[i], v)
                end
          a = out.boundaryField; b = u.boundaryField
          @simd for i in eachindex(a) 
                    @inbounds a[i] = $op(b[i], v)
                end
          out
      end 
    # memory-allocating versions
    @eval $op(u::ScalarField, v::Real) = $fname(u, v, similar(u))
end
# only multiplication is symmetric
*(v::Real, u::ScalarField) = u*v
mul!(v::Real, u::ScalarField, out::ScalarField) = mul!(u, v, out)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Overload operators on VectorFields ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# vector - vector
for (op, fname) in zip([:-, :+], [:sub!, :add!])
    # in-place versions
    @eval function $fname{D}(  u::VectorField{D}, 
                               v::VectorField{D}, 
                             out::VectorField{D})
              for d in 1:D 
                  $fname(u[d], v[d], out[d])
              end
              out
          end 
    # memory allocating versions
    @eval $op(u::VectorField, v::VectorField) = $fname(u, v, similar(u))
end

# vector - real
for (op, fname) in zip([:*, :/], [:mul!, :div!])
    @eval function $fname{D}(u::VectorField{D}, v::Real, out::VectorField{D})
             for d in 1:D
                 $fname(u[d], v, out[d])
             end
             out
          end 
    # memory allocating versions
    @eval $op(u::VectorField, v::Real) = $fname(u, v, similar(u)) 
end
# only multiplication is symmetric
*(v::Real, u::VectorField) = u*v
mul!(v::Real, u::VectorField, out::VectorField) = mul!(u, v, out)

""" Computes the quantity (u⋅∇)v - in-place.

    Notes
    -----
    In 3d this computes the vector field
    out[1] = u∂u/∂x + v*∂u/∂y + w*∂u/∂z
    out[2] = u∂v/∂x + v*∂v/∂y + w*∂v/∂z
    out[3] = u∂w/∂x + v*∂w/∂y + w*∂w/∂z
"""
function dotgrad!{D, T}( u::VectorField{D, T}, 
                        ∇v::TensorField{D, T}, 
                       out::VectorField{D, T}) 
    # set to zero the output
    fill!(out, zero(T))
    for di = 1:D
        for dj = 1:D
            muladd!(u[dj], ∇v[di, dj], out[di], out[di])
        end
    end
    out
end
# memory allocating version
dotgrad(u::VectorField, ∇v::TensorField) = dotgrad!(u, ∇v, similar(u))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Inner products, norms, and integrals ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Inner product between two vector fields """
@generated function inner{D, T, M}(u::VectorField{D, T, M}, 
                                   v::VectorField{D, T, M})
  # setup variables
    expr = quote
        I = zero(T)
        m = mesh(u)
        cvolumes_ = m.cvolumes
    end
    for d = 1:D
        ui, vi = symbol("u$d"), symbol("v$d")
        push!(expr.args, :($(ui) = u[$(d)].internalField))
        push!(expr.args, :($(vi) = v[$(d)].internalField))
    end
    # create loop part
    loop = :(@simd for i = 1:ncells(m) end)
    loopbody = loop.args[2].args[2].args
    push!(loopbody, :(@inbounds x = u1[i]*v1[i]))
    for d = 2:D
        ui, vi = symbol("u$d"), symbol("v$d")
        push!(loopbody, :(@inbounds x+= $(ui)[i] * $(vi)[i]))
    end
    push!(loopbody, :(@inbounds I += x*cvolumes_[i]))
    # now push the whole loop to expr and return 
    push!(expr.args, loop)
    push!(expr.args, :(return I))
    expr
end

""" Inner product between two scalar fields. This is 
    the integral of the product of the two fields. This
    is used for computing the integral of the product
    of two vorticity fields in 2D.
"""
function inner{D, T, M}(u::ScalarField{D, T, M}, v::ScalarField{D, T, M})
    I = zero(T)
    m = mesh(u)
    cvolumes_ = m.cvolumes
    ui_ = u.internalField
    vi_ = v.internalField
    @simd for i = 1:ncells(m)
        @inbounds I += ui_[i]*vi_[i]*cvolumes_[i]
    end
    I
end

""" L2 norm of vector or scalar field """
norm(u::Union{ScalarField, VectorField}) = sqrt(inner(u, u))

""" Compute partial derivative of `u` with respect to coordinate `dir`.
    
    Parameters
    ----------
      u : input scalar field
    out : scalar field containing ∂u/∂dir
    dir : either 1, 2, or 3

    Notes
    -----
    Gauss formula is used for computation of the average gradient in the 
    cell. A linear interpolation of the face value is used. This corresponds 
    to the 'Gauss linear' gradScheme and to the 'linear' interpolationScheme 
    in OpenFOAM.
"""
function der!{D, T}(u::ScalarField{D, T}, 
                    out::ScalarField{D, T}, 
                    dir::Integer)
    # for 2D simulations we cannot compute the gradient with respect 
    # to the third direction, it does not make sense to do it
    D == 2 && dir == 3 && error("Cannot compute partial derivative " *
                                "with respect to z for 2D field") 
    out.internalField[:] = zero(T)
    out.boundaryField[:] = zero(T)

    # need to get these fields out of the loop for better efficiency
    m = u.mesh
    αs = m.αs
    out_ = out.internalField
    ui_ = u.internalField
    ub_ = u.boundaryField
    fowners_ = m.fowners
    fneighs_ = m.fneighs
    fsvecs_ = m.fsvecs
    cvolumes_ = m.cvolumes

    # contributions only from non-empty boundary faces
    for patch in values(patches(m))
        if !(isempty(patch))
            for (faceID, ibnd) in patchfaces(m, patch)
                @inbounds out_[fowners_[faceID]] += (
                    getfield(fsvecs_[faceID], dir)*ub_[ibnd] )
            end
        end
    end

    # contributions from the internal faces
    @simd for faceID in internalfaces(m)
        @inbounds  begin 
            foi = fowners_[faceID]
            fni = fneighs_[faceID]
            fvalue = (1.0 - αs[faceID])*ui_[foi] + αs[faceID]*ui_[fni]
            value = getfield(fsvecs_[faceID], dir)*fvalue
            out_[foi] += value
            out_[fni] -= value
        end
    end

    # divide by the cell volume now, as for Gauss formula
    @simd for i = 1:ncells(m)
        out_[i] /= cvolumes_[i]
    end

    # FIXME: now we should fill the boundary field, 
    # either by interpolation or using the boundary conditions. This
    # requires some thoughts and programming.
    out
end

""" Compute gradient of `u` and write in `∇u`. """
function grad!{D}(u::VectorField{D}, ∇u::TensorField{D})
    for d = 1:D
        grad!(u[d], ∇u[d])
    end 
    ∇u
end
                        
""" Compute gradient of `u` and write in `∇u`. """
function grad!{D}(u::ScalarField{D}, ∇u::VectorField{D})
    for d = 1:D
        der!(u, ∇u[d], d)
    end 
    ∇u
end

""" Compute scalar vorticity of `u` and write in `ω`. Use `tmp` as storage. """
function curl!(u::VectorField{2}, ω::ScalarField{2}, tmp::ScalarField{2})
    der!(u[2], ω, 1)
    der!(u[1], tmp, 2)
    sub!(ω, tmp, ω)
end

# versions of the above which allocate the output
curl(u::VectorField) = curl!(u, ScalarField(mesh(u), ndims(u), eltype(u)), 
                                ScalarField(mesh(u), ndims(u), eltype(u)))
grad(u::ScalarField) = grad!(u, VectorField(mesh(u), ndims(u), eltype(u)))
grad(u::VectorField) = grad!(u, TensorField(mesh(u), ndims(u), eltype(u)))


# --- A few convenience functions ----
""" Compute average of the fields in `us`. """
function mean{T<:AbstractField}(us::AbstractVector{T})
    m = zero(us[1])
    for i = 1:length(us)
        add!(m, us[i], m)
    end
    mul!(m, 1.0/length(us), m)
end

""" Compute projections, (inner product), of each `VectorField` of `us` onto 
  each `VectorField` in `uis`. Returns a matrix containing:

  a[i, j] = inner(us[i], uis[j])

  An optional bias `m` can be given, for example the mean flow, which is 
  subtracted from each `VectorField` in `us` before projection.

"""
function projections{T<:VectorField}(us::AbstractVector{T}, 
                                     uis::AbstractVector{T};
                                     bias::Union{Bool, T}=false,
                                     verbose::Bool=true)
    M = length(us)
    N = length(uis)
    a = zeros(eltype(uis[1]), M, N)
    # need to subtract bias from us and write to tmp
    if bias != false
        tmp = similar(us[1])
        for i = 1:M
            tmp = sub!(us[i], bias, tmp)
            for j = 1:N
                a[i, j] = inner(tmp, uis[j])
            end
            verbose == true && print("\r Completed" *
              ": done $(round(100*i/M))%"); flush(STDOUT)
        end
    else 
        for i = 1:M
            for j = 1:N
                a[i, j] = inner(us[i], uis[j])
            end
            verbose == true && print("\r Completed" *
              ": done $(round(100*i/M))%"); flush(STDOUT)
        end
    end
    a
end

end