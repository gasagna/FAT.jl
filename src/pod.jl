# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module POD

import FAT.Fields: inner, mul!, add!

""" Snaphot POD algorithm 


    Parameters
    ----------
    u : a vector of objects that support an 'inner(u, v)' method
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
              " at: $(round(100.0*i/M, 1))%"); flush(STDOUT)
        for j = i:M
            val = convert(Float64, inner(u[i], u[j])) 
            C[i, j] = val/M
            C[j, i] = val/M
        end
    end

    # get eigen-decomposition
    λ, b = eig(C);
    λ = abs(flipdim(λ, 1))
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
                                 ": $(round(100.0*i/N, 1))%"); flush(STDOUT)
        # create a first zero field
        s = zero(u[1])
        # take a linear combination of the snapshots
        for j = 1:M
            mul!(u[j], b[j, i], tmp)
            add!(s, tmp, s)
            #s += u[j] * b[j, i]
        end
        push!(ui, mul!(s, 1.0/sqrt(M*λ[i]), s))
        ai[:, i] = b[:, i]*sqrt(M*λ[i])
    end
    verbose == true && println()
    ai, ui, λ
end

end