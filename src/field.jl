module Fields

import Base: call, *, /, +, -, norm, show, eltype, ndims, zero, mean

using FAT.Meshes
export AbstractField, ScalarField, VectorField, TensorField, inner, grad!, grad, 
       zeroField, curl!, zeroScalarField, zeroVectorField, zeroTensorField, 
       mesh, curl, integral, dotgrad!, dotgrad, projections


""" Abstract type for tensor fields

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

" Get the type of the field data "
eltype{D, T}(u::AbstractField{D, T}) = T

" Number of spatial dimensions "
ndims{D}(u::AbstractField{D}) = D

" Get the mesh of a field "
mesh(u::AbstractField) = u.mesh

# replace with generated functions, to replace the ntuple
# for scalar field the type parameters must be specified as D cannot be inferred
# from the input arguments
function zeroScalarField{T<:Real}(mesh::Mesh, D::Integer, dtype::Type{T}=Float64)
    D in [2, 3] || error("`D` must be either 2 or 3")
    ScalarField{D, dtype, typeof(mesh)}(zeros(dtype, ncells(mesh)), 
                                        zeros(dtype, nboundaryfaces(mesh)), 
                                        mesh)
end

zeroVectorField{T<:Real}(mesh::Mesh, D::Integer, dtype::Type{T}=Float64) =
    VectorField(ntuple(i->zeroScalarField(mesh, D, dtype), D), mesh)

zeroTensorField{T<:Real}(mesh::Mesh, D::Integer, dtype::Type{T}=Float64) =
    TensorField(ntuple(i->zeroVectorField(mesh, D, dtype), D), mesh)

# These are also useful
zero{D}(u::VectorField{D}) = zeroVectorField(mesh(u), D, eltype(u))
zero{D}(u::ScalarField{D}) = zeroScalarField(mesh(u), D, eltype(u))

# overload arithmetic operators on scalar and vector fields. These
# should not be used in performance critical code, because lots of 
# temporaries are created. 
*{D}(u::ScalarField{D}, a::Real) = ScalarField{D, eltype(u), typeof(mesh(u))}(scale(u.internalField, a), scale(u.boundaryField, a), mesh(u))   
*{D}(u::VectorField{D}, a::Real) = VectorField(ntuple(i -> u.scalars[i]*a, D), mesh(u))
*(a::Real, u::ScalarField) = u*a   
*(a::Real, u::VectorField) = u*a
/(u::ScalarField, a::Real) = u*inv(a)
/{D}(u::VectorField{D}, a::Real) = VectorField(ntuple(i -> u.scalars[i]*inv(a), D), mesh(u))

# operations on Fields of equal rank. We assume they have 
# same mesh and same amount of data. It would be strange 
# if they did not, but we would have some error here
for op in [:-, :+]
    @eval $op{D, T, M}(u::ScalarField{D, T, M}, v::ScalarField{D, T, M}) = 
        ScalarField{D, T, M}($op(u.internalField, v.internalField), $op(u.boundaryField, v.boundaryField), mesh(u))
    @eval $(op){D}(u::VectorField{D}, v::VectorField{D}) = 
        VectorField(ntuple(i -> $op(u.scalars[i], v.scalars[i]), D), mesh(u))
end 

# fast in-place dot product of two vector fields
function mul!{D}(u::VectorField{D}, v::VectorField{D}, out::ScalarField{D})
    # set out to 0 before doing anything else
    out.internalField[:] = zero(eltype(out))
    out.boundaryField[:] = zero(eltype(out))
    for d = 1:D
        # muladd!(u, v, w, out) -> out = u*v + w for ScalarFields
        muladd!(u.scalars[d], v.scalars[d], out, out)
    end
    out
end
# memory-allocating version
*{D}(u::VectorField{D}, v::VectorField{D}) = 
    mul!(u, v, zeroScalarField(mesh(u), D, eltype(u)))

# These are fast, in-place operations on fields, and should be used in performance critical code
# operations with scalars
for (op, fname) in zip([:-, :+, :*, :/], [:sub!, :add!, :mul!, :div!])
    @eval function $fname(u::ScalarField, v::ScalarField, out::ScalarField)
             @inbounds begin
                 a = out.internalField 
                 b = u.internalField         
                 c = v.internalField  
                 @simd for i in eachindex(u.internalField)
                           a[i] = $op(b[i], c[i])
                       end
                 a = out.boundaryField
                 b = u.boundaryField
                 c = v.boundaryField
                 @simd for i in eachindex(u.boundaryField) 
                           a[i] = $op(b[i], c[i])
                       end
             end
             out
          end 
    @eval function $fname(u::ScalarField, v::Real, out::ScalarField)
             @inbounds begin 
                 a = out.internalField 
                 b = u.internalField         
                 @simd for i in eachindex(u.internalField)
                           a[i] = $op(b[i], v)
                       end
                 a = out.boundaryField
                 b = u.boundaryField
                 @simd for i in eachindex(u.boundaryField)
                           a[i] = $op(b[i], v)
                       end
             end
             out
          end 
    @eval $fname(v::Real, u::ScalarField, out::ScalarField) = $fname(u, v, out)
end   
# memory-allocating versions
*{D}(u::ScalarField{D}, v::ScalarField{D}) = mul!(u, v, zero(u))

# compute out = u*v + w efficiently
function muladd!(u::ScalarField, v::ScalarField, w::ScalarField, out::ScalarField)
    @inbounds @simd for i in eachindex(u.internalField)
        out.internalField[i] = u.internalField[i] * v.internalField[i] + w.internalField[i]
    end
    @inbounds @simd for i in eachindex(u.boundaryField)
        out.boundaryField[i] = u.boundaryField[i] * v.boundaryField[i] + w.boundaryField[i]
    end
    out
end   

# operations between vectors
for fname in [:sub!, :add!]
    @eval function $fname{D}(u::VectorField{D}, v::VectorField{D}, out::VectorField{D})
              for d in 1:D 
                  $fname(u.scalars[d], v.scalars[d], out.scalars[d])
              end
              out
          end 
end

# operations between vector fields and a scalar
for fname in [:sub!, :add!, :mul!, :div!]
    @eval function $fname{D}(u::VectorField{D}, v::Real, out::VectorField{D})
             for d in 1:D
                 $fname(u.scalars[d], v, out.scalars[d])
             end
             out
          end 
    @eval $fname{D}(v::Real, u::VectorField{D}, out::VectorField{D}) = $fname(u, v, out)
end

""" Compute the term (u⋅∇)v and write it in `out`.

    Notes
    -----
    In 3d this computes the vector field
    out[1] = u∂u/∂x + v*∂u/∂y + w*∂u/∂z
    out[2] = u∂v/∂x + v*∂v/∂y + w*∂v/∂z
    out[3] = u∂w/∂x + v*∂w/∂y + w*∂w/∂z

    This is an in-place version. A version allocating the output also exists.
"""
function dotgrad!{D}(u::VectorField{D}, ∇v::TensorField{D}, out::VectorField{D}) 
    for di = 1:D
        # set to zero the output
        out.scalars[di].internalField[:] = zero(eltype(u))
        for dj = 1:D
            muladd!(u.scalars[dj], ∇v.vectors[di].scalars[dj], out.scalars[di], out.scalars[di])
        end
    end
    out
end
dotgrad{D}(u::VectorField{D}, ∇v::TensorField{D}) = dotgrad!(u, ∇v, zero(u))


""" Inner product between two vector fields """
function inner{D, T, M}(u::VectorField{D, T, M}, v::VectorField{D, T, M})
    I = zero(T)
    m = mesh(u)
    c = m.cells
    for d = 1:D
        ud = u.scalars[d].internalField
        vd = v.scalars[d].internalField
        @simd for i = 1:ncells(m)
            @inbounds I += ud[i]*vd[i]*volume(c[i])
        end
    end
    I
end

""" Inner product between two scalar fields. This is the integral 
    of the product of the two fields. 
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

""" Compute integral of scalar field `u`. """
function integral{D, T}(u::ScalarField{D, T}) 
    I = zero(T)
    @inbounds begin 
        for (i, cell) in enumerate(cells(mesh(u)))
            I += u.internalField[i]*volume(cell)
        end
    end
    I
end

""" L2 norm of vector field """
norm(u::VectorField) = sqrt(inner(u, u))

""" L2 norm of scalar field """
norm(u::ScalarField) = sqrt(integral(mul!(u, u, u)))

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













