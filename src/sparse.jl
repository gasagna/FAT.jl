# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

module Sparse

export monoexps, 
       nmonomials, 
       regrmat, 
       regrmat!

"""
    Calculates the exponents of the state variables in each 
    monomial term for a quadratic polynomial in `N` variables.

    Arguments
    ---------
    N::Integer  : the number of state variables

    Returns
    -------
    out::Vector{NTuple{N, Int}} : a vector of `1 + N + N*(N+1)/2` 
        tuples with the exponents of the i-the monomial.
"""
function monoexps(N::Integer)
    # allocate vector of tuples
    exps = Vector{NTuple{N, Int}}()
    sizehint!(exps, nmonomials(N))

    # constant term
    push!(exps, ntuple(j->0, N))

    # linear term
    for j in 1:N
        push!(exps, ntuple(i->i==j ? 1 : 0, N))
    end
    
    # quadratic term
    for j in 1:N, k in j:N
        push!(exps, ntuple(i->i==j && i==k ? 2 :
                              i==k ? 1 :
                              i==j ? 1 : 0, N))
    end
    exps
end

" Number of monomials in a quadratic multivariate polynomial in `N` variables. "
nmonomials(N::Integer) = 1 + N + N*(N+1)>>1

""" Compute data matrix `R` of regressors from time histories of state vector
    in the matrix `A`. The matrix `A` is assumed to be organised such that the
    `i`-th row is the time history of state variable `i`. This arrangement
    improves the memory access pattern 

    Arguments
    ---------
    A::Matrix : the data matrix

    Returns
    -------
    R::Matrix : the regressors matrix
"""
regrmat{T<:Real}(A::Matrix{T}) = regrmat!(A, ones(T, nmonomials(size(A, 1)), size(A, 2)))

# inplace version
function regrmat!{T<:Real}(A::AbstractMatrix{T}, R::AbstractMatrix{T})
    nrowsA, ncolsA = size(A)
    nrowsR, ncolsR = size(R)
    ncolsR == ncolsR || error("incompatible number of columns of `A` and `R`")
    nrowsR == nmonomials(nrowsA) || error("incompatible number of rows of `A` and `R`")
    nrowsA > ncolsA && error("matrix has more columns than rows.. r u sure?")
    for col = 1:ncolsA
        regrvec!(slice(A, :, col), slice(R, :, col))
    end
    R
end

"""
    Compute regressors vector from state vector `a = [a_1 ... a_N]`.
    
    Arguments
    ---------
    a::Vector : `N`-dimensional state vector
    r::Vector : `nmonomials(N)`-dimensional vector of regressors
"""
function regrvec!{T<:Real}(a::AbstractVector{T}, r::AbstractVector{T})
    N = length(a)
    @inbounds begin
        # constant term
        r[1] = one(T)
        # linear term
        for i = 1:N
            r[i+1] = a[i]
        end
        # quadratic terms
        curr = N+2
        for i in 1:N, j in i:N
            r[curr] = a[i]*a[j]
            curr += 1
        end
    end
    r
end
end