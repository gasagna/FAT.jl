# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Fields

import Base: call, 
             *, 
             /, 
             +,
             -,
             norm,
             show,
             eltype,
             ndims,
             zero,
             mean,
             fill!,
             similar

using FAT.Meshes
export AbstractField, ScalarField, VectorField, TensorField, inner, grad!
export curl!, grad, mesh, curl, integral, dotgrad!, dotgrad
export projections


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
    show(io, mesh(u); space="     ")
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


" Get the type of the field data "
eltype{D, T}(u::AbstractField{D, T}) = T

" Number of spatial dimensions "
ndims{D}(u::AbstractField{D}) = D

" Get the mesh of a field "
mesh(u::AbstractField) = u.mesh


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
        fill!(u.scalars[d], zero(T))
    end
    u
end
zero{D, T}(u::VectorField{D, T}) = fill!(VectorField(mesh(u), D, T), zero(T))

function fill!{D, T}(u::TensorField{D, T}, val::Real)
    for d = 1:D
        fill!(u.vectors[d], zero(T))
    end
    u
end
zero{D, T}(u::TensorField{D, T}) = fill!(TensorField(mesh(u), D, T), zero(T))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Overload operators on ScalarFields ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note we never check if dimensions and mesh of the operand match
for (op, fname) in zip([:-, :+, :*], [:sub!, :add!, :mul!])
    # fast in-place version
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
    c = v.internalField;   d = v.internalField  
    @simd for i in eachindex(a)
        @inbounds a[i] = b[i]*c[i] + d[i]
    end
    a = out.boundaryField; b = u.boundaryField
    c = v.boundaryField;   d = v.boundaryField  
    @simd for i in eachindex(a)
        @inbounds a[i] = b[i]*c[i] + d[i]
    end
    out
end   

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Overload operators on VectorFields ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (op, fname) in zip([:-, :+], [:sub!, :add!])
    # in-place versions
    @eval function $fname{D}(  u::VectorField{D}, 
                               v::VectorField{D}, 
                             out::VectorField{D})
              for d in 1:D 
                  $fname(u.scalars[d], v.scalars[d], out.scalars[d])
              end
              out
          end 
    # memory allocating versions
    @eval $op(u::VectorField, v::VectorField) = $fname(u, v, similar(u))
end

""" Compute the term (u⋅∇)v and write it in `out`.

    Notes
    -----
    In 3d this computes the vector field
    out[1] = u∂u/∂x + v*∂u/∂y + w*∂u/∂z
    out[2] = u∂v/∂x + v*∂v/∂y + w*∂v/∂z
    out[3] = u∂w/∂x + v*∂w/∂y + w*∂w/∂z

    This is an in-place version.
"""
function dotgrad!{D, T}( u::VectorField{D, T}, 
                        ∇v::TensorField{D, T}, 
                       out::VectorField{D, T}) 
    # set to zero the output
    fill!(out, zero(T))
    for di = 1:D
        for dj = 1:D
            muladd!(u.scalars[dj], 
                    ∇v.vectors[di].scalars[dj], 
                    out.scalars[di], out.scalars[di])
        end
    end
    out
end
# memory allocating version
dotgrad{D}(u::VectorField{D}, ∇v::TensorField{D}) = dotgrad!(u, ∇v, similar(u))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Inner products, norms, and integrals ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Inner product between two vector fields """
@generated function inner{D, T, M}(u::VectorField{D, T, M}, v::VectorField{D, T, M})
  # setup variables
    expr = quote
        I = zero(T)
        m = mesh(u)
        c = m.cells
    end
    for d = 1:D
        ui, vi = symbol("u$d"), symbol("v$d")
        push!(expr.args, :($(ui) = u.scalars[$(d)].internalField))
        push!(expr.args, :($(vi) = v.scalars[$(d)].internalField))
    end
    # create loop part
    loop = :(@simd for i = 1:ncells(m) end)
    loopbody = loop.args[2].args[2].args
    push!(loopbody, :(@inbounds x = u1[i]*v1[i]))
    for d = 2:D
        ui, vi = symbol("u$d"), symbol("v$d")
        push!(loopbody, :(@inbounds x+= $(ui)[i] * $(vi)[i]))
    end
    push!(loopbody, :(@inbounds I += x*volume(c[i])))
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
    c = m.cells
    ud = u.internalField
    vd = v.internalField
    @simd for i = 1:ncells(m)
        @inbounds I += ud[i]*vd[i]*volume(c[i])
    end
    I
end

""" L2 norm of vector field """
norm(u::VectorField) = sqrt(inner(u, u))

""" L2 norm of scalar field """
norm(u::ScalarField) = sqrt(inner(u, u))

""" Compute partial derivative of `u` with respect to coordinate `dir`.
    
    Parameters
    ----------
      u : input scalar field
    out : scalar field containing ∂u/∂dir
    dir : either 1, 2, or 3, for x, y, or z

    Notes
    -----
    Gauss formula is used for computation of the average gradient in the cell. A linear
    interpolation of the face value is used. This corresponds to the 'Gauss linear' 
    gradScheme and to the 'linear' interpolationScheme in OpenFOAM.
"""
function der!{D, T}(u::ScalarField{D, T}, out::ScalarField{D, T}, dir::Union{Integer, Symbol})
    # for 2D simulation we cannot compute the gradient with respect to the third direction
    D == 2 && (dir == 3 || dir == :z) && error("Cannot compute partial "*
                                                "derivative with respect to z for 2D field") 
    @inbounds begin
        out.internalField[:] = zero(T)
        out.boundaryField[:] = zero(T)
        # FIXME we should only loop over the non-empty patches. This
        # does not do any harm now, because empty face have 'usually' opposite 
        # surface vectors that cancel each other. In case of axisymmetric problems
        # this should be fixed. However, this is inefficient.
        for (i, face) in enumerate(boundaryfaces(mesh(u)))
            # println(i, " ", svec(face, Val{dir}), " ", u.boundaryField[i])
            out.internalField[ownerID(face)] += svec(face, Val{dir})*u.boundaryField[i]
        end
        for (α, face) in zip(mesh(u).αs, internalfaces(mesh(u)))
            face_value = (one(α)-α)*u.internalField[ownerID(face)] + α*u.internalField[neighbourID(face)]
            out.internalField[ownerID(face)]     += svec(face, Val{dir})*face_value
            out.internalField[neighbourID(face)] -= svec(face, Val{dir})*face_value
        end
        for (i, cell) in enumerate(cells(mesh(u)))
            out.internalField[i] /= volume(cell)
        end
        # FIXME: now we should fill the boundary field, 
        # either by interpolation or using the boundary conditions
    end
    out
end

""" Compute gradient of `u` and write in `∇u`. """
function grad!{D}(u::VectorField{D}, ∇u::TensorField{D})
    for d = 1:D
        grad!(u.scalars[d], ∇u.vectors[d])
    end 
    ∇u
end
                        
""" Compute gradient of `u` and write in `∇u`. """
function grad!{D}(u::ScalarField, ∇u::VectorField{D})
    for d = 1:D
        der!(u, ∇u.scalars[d], d)
    end 
    ∇u
end

""" Compute scalar vorticity of `u` and write in `ω`. Use `tmp` as temporary storage. """
function curl!(u::VectorField{2}, ω::ScalarField{2}, tmp::ScalarField{2})
    der!(u.scalars[2], ω, 1)
    der!(u.scalars[1], tmp, 2)
    sub!(ω, tmp, ω)
end

# versions of the above which allocate the output
curl(u::VectorField) = curl!(u, zeroScalarField(mesh(u), ndims(u), eltype(u)), 
                                zeroScalarField(mesh(u), ndims(u), eltype(u)))
grad(u::ScalarField) = grad!(u, zeroVectorField(mesh(u), ndims(u), eltype(u)))
grad(u::VectorField) = grad!(u, zeroTensorField(mesh(u), ndims(u), eltype(u)))


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
                                     bias::Union{Bool, T}=false)
    M = length(us)
    N = length(uis)
    a = zeros(eltype(uis[1]), M, N)
    # need to subtract bias from us and write to tmp
    if bias != false
        tmp = zero(us[1])
        for i = 1:M
            tmp = sub!(us[i], bias, tmp)
            for j = 1:N
                a[i, j] = inner(tmp, uis[j])
            end
        end
    else 
        for j = 1:N, i = 1:M
            a[i, j] = inner(us[i], uis[j])
        end
    end
    a
end

end













