using Base.Test
using Benchmarks

include("sparse.jl")

# test number of exponents
@test nmonomials(0) == 1
@test nmonomials(1) == 3
@test nmonomials(2) == 6

for i = 5:5:200
    s = @sprintf "%d" i 
    println(lpad(s, 4), " : ", nmonomials(i))
end


# test exponents
exps = monoexps(4)
@test exps[1]  == (0, 0, 0, 0)
@test exps[2]  == (1, 0, 0, 0)
@test exps[3]  == (0, 1, 0, 0)
@test exps[4]  == (0, 0, 1, 0)
@test exps[5]  == (0, 0, 0, 1)
@test exps[6]  == (2, 0, 0, 0)
@test exps[7]  == (1, 1, 0, 0)
@test exps[8]  == (1, 0, 1, 0)
@test exps[9]  == (1, 0, 0, 1)
@test exps[10] == (0, 2, 0, 0)
@test exps[11] == (0, 1, 1, 0)
@test exps[12] == (0, 1, 0, 1)
@test exps[13] == (0, 0, 2, 0)
@test exps[14] == (0, 0, 1, 1)
@test exps[15] == (0, 0, 0, 2)
@test length(exps) == 15


# test generation of regressors vector 
a = [1, 2]
r = zeros(eltype(a), nmonomials(length(a)))
r = regrvec!(a, r)
@test r[1] == 1
@test r[2] == 1
@test r[3] == 2
@test r[4] == 1
@test r[5] == 2
@test r[6] == 4

# test generation of regressor matrix
A = repmat(a, 1, 3)
R = regrmat(A)
@test R == [1 1 1;
            1 1 1;
            2 2 2;
            1 1 1; 
            2 2 2;
            4 4 4]