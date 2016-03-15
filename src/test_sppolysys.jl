include("sppolysys.jl")
using FAT.Sparse



N = 2
M = nmonomials(N)
println(N, " ", M)
c = randn(M)

# println(macroexpand(quote @sparsePolySys N c... end))

function fgen(N, c) 
	f(x) = @sparsePolySys(N, c...)
end

x = randn(N)

println(f(x))

println(@code_llvm f(x))
