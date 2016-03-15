# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #

# macro sparsePolySys(N, coeffs...)
#     def = Expr(:block)
#     for i = 1:N
#         xi = symbol("x_$i")
#         push!(def.args, :(@inbounds $xi = x[$i]))
#     end
#     # constant term
#     exec = esc(coeffs[1])
#     curr = 2
#     # linear term
#     for j = 1:N
#         Lij = esc(coeffs[curr])
#         xj = symbol("x_$j")
#         exec = :($exec + $Lij*$xj)
#         curr += 1
#     end
#     # quadratic term
#     for j = 1:N, k = j:N
#         Qijk = esc(coeffs[curr])
#         xj = symbol("x_$j")
#         xk = symbol("x_$k")
#         exec = :($exec + $Qijk*$xj*$xk)
#         curr += 1
#     end
#     Expr(:block, def, exec)
# end




type SparsePolynomialSystem{T, N}
  # Vector of coefficients
  θ::Vector{T}
  # non zero entries of constant, linear and quadratic parts
  nzi::Vector{NTuple{1, T}}
  nzij::Vector{NTuple{2, T}}
  nzijk::Vector{NTuple{3, T}}
  istart::Vector{Int}
end

for (i, sym,  nzvec) in zip([1, 2, 3], [:C, :L, :N], [:nzi, :nzij, :nzijk])
  @eval begin
      immutable $sym
          sys::SparsePolynomialSystem
      end
      length(s::$sym) = length(s.sys.$nzvec)
      start(s::$sym) = 1
      next(s::$sym, state) = (s.sys.$nzvec[state], s.sys.θ[s.sys.istart[$i] + state]), state+1
      done(s::$sym, state) = state == length(s.sys.$nzvec)
  end
end

function f!{T}(sys::SparsePolynomialSystem{T}, t::T, x::Vector{T}, out::Vector{T})
  # set to zero first
  fill!(out, zero(T))
  for (i, Ci) in C(sys) 
      out[i] += Ci
  end
  for ((i, j), Lij) in L(sys)
      out[i] += Lij*x[j]
  end
  for ((i, j, k), Nijk) in N(sys)
      out[i] += Nijk*x[j]*x[k]
  end
end