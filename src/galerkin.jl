# ------------------------------------------------------------------- #
# Copyright 2015-2016, Davide Lasagna, AFM, University of Southampton #
# ------------------------------------------------------------------- #
module Galerkin

export GalerkinModel, tofile, fromfile

import JLD: save,
			load

import Base: call

import FAT.Fields: inner,
				   VectorField,
				   curl,
				   grad,
				   dotgrad!

type GalerkinModel
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

# save model to a .jld file
tofile(filename::AbstractString, sys::GalerkinModel) = 
	save(filename, "sys", sys)

# load model from a .jld file
fromfile(filename::AbstractString) = load(filename, "sys")

""" Build a Galerkin model from the basis functions `us` and base 
	flow `u0`, for Reynolds number `Re`.

	Parameters
	----------
	  u0 : the base flow
	  us : the basis functions
	  Re : Reynolds number
	symm : enforce symmetries in the nonlinear term
"""
function GalerkinModel{T<:VectorField}(u0::T, 
									   uis::Vector{T}, 
									   Re::Real; symm::Bool=true)
	# Precompute gradients and vorticity fields
	ω0 = curl(u0)
	∇u0 = grad(u0)
	ωis = map(curl, uis)
	∇uis = map(grad, uis)

	# build system
	N = length(uis)
	sys = GalerkinModel(N)

	# temporary vector field
	u∇u = zero(u0)

	# constant term
	for i in 1:N
		sys.c[i]  = ( - inner(ωis[i], ω0)/Re
		              - inner(uis[i], dotgrad!(u0, ∇u0, u∇u)) )
	end

	# linear term 
	for i = 1:N, j = 1:N
		sys.L[i, j] = ( - inner(ωis[i], ωis[j])/Re 
					    - inner(uis[i], dotgrad!(u0, ∇uis[j], u∇u))
					    - inner(uis[i], dotgrad!(uis[j], ∇u0, u∇u)) )
	end

	# quadratic term. This could save some symmetries
	for k = 1:N, j = 1:N, i = 1:N
		if i != k 
			sys.Q[i, j, k] = - inner(uis[i], dotgrad!(uis[j], ∇uis[k], u∇u))
		end
	end

	if symm == true
		# manually enforce symmetries in the nonlinear term
		for i = 1:N, j = 1:N, k = 1:N
			a = sys.Q[i, j, k]
			b = sys.Q[k, j, i]
			sys.Q[i, j, k] = abs((a-b)/2)*sign(a)
			sys.Q[k, j, i] = abs((a-b)/2)*sign(b)
		end
	end
	sys
end


# Galerkin models are callable objects.
function call(sys::GalerkinModel, xdot::AbstractVector, x::AbstractVector)
	N = length(x)
	@inbounds begin 
		for i = 1:N
			xdot[i] = sys.c[i] 
		end
		for j = 1:N, i = 1:N
			xdot[i] += sys.L[i, j]*x[j] 
		end
		for k = 1:N, j = 1:N, i = 1:N
			xdot[i] += sys.Q[i, j, k]*x[j]*x[k]
		end
	end
	xdot
end

end