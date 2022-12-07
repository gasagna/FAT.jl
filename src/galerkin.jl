# ------------------------------------------------------------------- #
# Copyright 2015-2019, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Galerkin

export GalerkinModel

import FAT.Fields: ScalarField,
                   VectorField,
                   TensorField,
                   curl,
                   grad,
                   dotgrad!,
                   mesh

struct GalerkinModel
    c::Array{Float64, 1}
    L::Array{Float64, 2}
    Q::Array{Float64, 3}
    # Create an empty type, with zeroes, that will be filled later
    function GalerkinModel(N::Integer)
        c = zeros(Float64, N)
        L = zeros(Float64, N, N)
        Q = zeros(Float64, N, N, N)
        new(c, L, Q)
    end
end

""" Build a Galerkin model from the basis functions `us` and base
    flow `u0`, for Reynolds number `Re`.

    Parameters
    ----------
      u0 : the base flow
      us : the basis functions
      Re : Reynolds number
    symm : enforce symmetries in the nonlinear term
"""
function GalerkinModel(u0::VF,
                      uis::Vector{VF},
                       Re::Real;
                     symm::Bool=true,
                  verbose::Bool=true) where {VF<:VectorField}
    # build system
    N = length(uis)
    sys = GalerkinModel(N)

    # Precompute gradients and vorticity fields
    ω0 = curl(u0)
    ∇u0 = grad(u0)
    ωis  = []
    ∇uis = []
    for i = 1:N
        gr = grad(uis[i])
        push!(∇uis, gr)
        push!(ωis, gr[2, 1] .- gr[1, 2])
        verbose == true && print("\r Calculation of velocity gradient fields" *
              ": done $i over $N"); flush(stdout)
    end
    verbose == true && println()

    # temporary vector field
    u∇u = zero(u0)

    # constant term
    for i in 1:N
        sys.c[i] =  - (ωis[i] ⋅ ω0)/Re - uis[i] ⋅ dotgrad!(u0, ∇u0, u∇u)
        verbose == true && print("\r Constant term" *
              ": done $(round(100*i/N, 1))%"); flush(stdout)
    end
    verbose == true && println()


    # linear term
    for i = 1:N
        for j = 1:N
            sys.L[i, j] = ( - (ωis[i] ⋅ ωis[j])/Re
                            -  uis[i] ⋅ dotgrad!(u0, ∇uis[j], u∇u)
                            -  uis[i] ⋅ dotgrad!(uis[j], ∇u0, u∇u))
        end
        verbose == true && print("\r Linear term" *
              ": done $(round(100*i/N, 1))%"); flush(stdout)
    end
    verbose == true && println()


    # quadratic term. This could save some symmetries
    for k = 1:N
        for j = 1:N, i = 1:N
            if i != k
                sys.Q[i, j, k] = - uis[i] ⋅ dotgrad!(uis[j], ∇uis[k], u∇u)
                println(i, " ", j, " ", k)
            end
        end
        verbose == true && print("\r Nonlinear term" *
              ": done $(round(100*k/N, 1))%"); flush(stdout)
    end
    verbose == true && println()


    if symm == true
        # manually enforce symmetries in the nonlinear term
        for i = 1:N, j = 1:N, k = 1:N
            a = sys.Q[i, j, k]
            b = sys.Q[k, j, i]
            sys.Q[i, j, k] = abs((a-b)/2)*sign(a)
            sys.Q[k, j, i] = abs((a-b)/2)*sign(b)
        end
    end
    return sys
end


# Galerkin models are callable objects.
function (sys::GalerkinModel)(t::Real, x::AbstractVector, xdot::AbstractVector)
    N = length(x)
    @inbounds begin
        @simd for i = 1:N
            xdot[i] = sys.c[i]
        end
        for j = 1:N
            @simd for i = 1:N
                xdot[i] += sys.L[i, j]*x[j]
            end
        end
        for k = 1:N, j = 1:N
            @simd for i = 1:N
                xdot[i] += sys.Q[i, j, k]*x[j]*x[k]
            end
        end
    end
    return xdot
end

end