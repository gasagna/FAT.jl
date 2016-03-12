using Base.Test
using FAT.Galerkin

# create a random model
srand(1)
N = 10
sys = GalerkinModel(N)
sys.c = rand(N)
sys.L = rand(N, N)
sys.Q = rand(N, N, N)

# save it
tofile("here.jld", sys)

# load it
sys2 = fromfile("here.jld")

@test sys2.c == sys.c
@test sys2.L == sys.L
@test sys2.Q == sys.Q

# clean up
rm("here.jld")

