# ------------------------------------------------------------------- #
# Copyright 2015-2019, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

""" Computes the quantity (u⋅∇)v - in-place.

    Notes
    -----
    In 3d this computes the vector field
    out[1] = u∂u/∂x + v*∂u/∂y + w*∂u/∂z
    out[2] = u∂v/∂x + v*∂v/∂y + w*∂v/∂z
    out[3] = u∂w/∂x + v*∂w/∂y + w*∂w/∂z
"""
function dotgrad!(u::VectorField{2}, ∇v::TensorField{2}, out::VectorField{2})
    out[1] .= u[1] .* ∇v[1, 1] .+ u[2] .* ∇v[1, 2]
    out[2] .= u[1] .* ∇v[2, 1] .+ u[2] .* ∇v[2, 2]
    return out
end

# memory allocating version
dotgrad(u::VectorField, ∇v::TensorField) = dotgrad!(u, ∇v, similar(u))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ Inner products, norms, and integrals ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""" Inner product between two vector fields """
LinearAlgebra.dot(u::VF, v::VF) where {D, VF<:VectorField{D}} =
    sum(u.scalars[d] ⋅ v.scalars[d] for d = 1:D)

""" Inner product between two scalar fields. This is 
    the integral of the product of the two fields. This
    is used for computing the integral of the product
    of two vorticity fields in 2D.
"""
function LinearAlgebra.dot(u::SF, v::SF) where {SF<:ScalarField}
    I = zero(eltype(u))
    m = mesh(u)
    cvolumes_ = m.cvolumes
    ui_ = u.internalField
    vi_ = v.internalField
    @simd for i = 1:ncells(m)
        @inbounds I += ui_[i]*vi_[i]*cvolumes_[i]
    end
    return I
end

""" L2 norm of vector or scalar field """
LinearAlgebra.norm(u::Union{ScalarField, VectorField}) = sqrt(u ⋅ u)

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
function der!(  u::ScalarField{D, T}, 
              out::ScalarField{D, T}, 
               dir::Integer) where {D, T}
    # for 2D simulations we cannot compute the gradient with respect 
    # to the third direction, it does not make sense to do it
    D == 2 && dir == 3 && error("Cannot compute partial derivative " *
                                "with respect to z for 2D field") 
    out.internalField .= zero(T)
    out.boundaryField .= zero(T)

    # aliases
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
    for (patchname, patch) in patches(m)
        if !(isempty(patch))
            for (faceID, ibnd) in faceiterator(m, patchname)
                @inbounds out_[fowners_[faceID]] += (
                    getfield(fsvecs_[faceID], dir)*ub_[ibnd] )
            end
        end
    end

    # contributions from the internal faces
    @simd for faceID in 1:ninternalfaces(m)
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
    out_ ./= cvolumes_

    # FIXME: now we should fill the boundary field of the derivative 
    # either by interpolation or using the boundary conditions. This
    # requires some thoughts and programming.
    return out
end

""" Compute gradient of `u` and write in `∇u`. """
function grad!(u::VectorField{D}, ∇u::TensorField{D}) where {D}
    for d = 1:D
        grad!(u[d], ∇u[d])
    end 
    return ∇u
end
                        
""" Compute gradient of `u` and write in `∇u`. """
function grad!(u::ScalarField{D}, ∇u::VectorField{D}) where {D}
    for d = 1:D
        der!(u, ∇u[d], d)
    end 
    return ∇u
end

""" Compute scalar vorticity of `u` and write in `ω`. Use `tmp` as storage. """
function curl!(u::VectorField{2}, ω::ScalarField{2}, tmp::ScalarField{2})
    ω .= der!(u[2], ω, 1) .- der!(u[1], tmp, 2)
    return ω
end

# versions of the above which allocate the output
curl(u::VectorField) = curl!(u, ScalarField(mesh(u), ndims(u), eltype(u)), 
                                ScalarField(mesh(u), ndims(u), eltype(u)))
grad(u::ScalarField) = grad!(u, VectorField(mesh(u), ndims(u), eltype(u)))
grad(u::VectorField) = grad!(u, TensorField(mesh(u), ndims(u), eltype(u)))


# --- A few convenience functions ----
""" Compute average of the fields in `us`. """
function Statistics.mean(us::AbstractVector{<:AbstractField})
    m = copy(us[1])
    for i = 2:length(us)
        m .+= us[i]
    end
    m ./= length(us)
    return m
end

""" Compute projections, (dot product), of each `VectorField` of `us` onto 
  each `VectorField` in `uis`. Returns a matrix containing:

  a[i, j] = dot(us[i], uis[j])

  An optional bias `m` can be given, for example the mean flow, which is 
  subtracted from each `VectorField` in `us` before projection.

"""
function projections( us::AbstractVector{T}, 
                     uis::AbstractVector{T};
                    bias::Union{Bool, T}=false,
                 verbose::Bool=true) where {T<:VectorField}
    M = length(us)
    N = length(uis)
    a = zeros(eltype(uis[1]), M, N)
    # need to subtract bias from us and write to tmp
    if bias != false
        tmp = similar(us[1])
        for i = 1:M
            tmp .= us[i] .- bias
            for j = 1:N
                a[i, j] = dot(tmp, uis[j])
            end
            verbose == true && print("\r Completed" *
              ": done $(round(100*i/M))%"); flush(stdout)
        end
    else 
        for i = 1:M
            for j = 1:N
                a[i, j] = dot(us[i], uis[j])
            end
            verbose == true && print("\r Completed" *
              ": done $(round(100*i/M))%"); flush(stdout)
        end
    end
    return a
end
