# ------------------------------------------------------------------- #
# Copyright 2015-2022, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module POD

""" Snaphot POD algorithm 


    Parameters
    ----------
    u : a vector of objects that support an 'dot(u, v)' method
    N : the desired number of POD modes
"""
function snapshotPOD(u::AbstractVector, N::Integer; verbose::Bool=true)
    # number of snapshots
    M = length(u)

    # correlation matrix
    C = Array(Float64, M, M)

    # compute it
    for i = 1:M
        verbose == true && print("\r Correlation matrix completed" *
              " at: $(round(100.0*i/M, 1))%"); flush(stdout)
        for j = i:M
            val = u[i] ⋅ u[j]
            C[i, j] = val/M
            C[j, i] = val/M
        end
    end

    # get eigen-decomposition
    λ, b = eig(LinearAlgebra.Symmetric(C));
    λ = flipdim(λ, 1)
    b = flipdim(b, 2)

    # initialise
    ui = eltype(u)[]
    ai = zeros(M, N)

    # temporary vector field
    tmp = zero(u[1])

    verbose == true && println()

    # for each mode we need
    for i = 1:N
        verbose == true && print("\r Construction of modes" *
                                 ": $(round(100.0*i/N; digits=1))%"); flush(stdout)
        # first create an empty field
        s = zero(u[1])

        # take a linear combination of the snapshots
        for j = 1:M
            s .+= u[j] .* b[j, i]
        end
        
        # normalise
        s .*= 1/sqrt(M*λ[i])
        @views ai[:, i] .= b[:, i]*sqrt(M*λ[i])
        
        # store
        push!(ui, s)
    end
    verbose == true && println()
    return ai, ui, λ
end

end