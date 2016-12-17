using Base.Test

using FAT.Fields
using FAT.Meshes
using FAT.Simulations
using FAT.Utils

casedirname = "ldc_test_binary"

# use a test run
sim = SimulationData(casedirname)

# check times
@test times(sim) == collect(0:10)

# check we return a mesh
@test typeof(mesh(sim)) == Mesh{Float64}

# check casedirname
@test casedir(sim) == casedirname

# # check loading snapshots with the load_snapshot function, which should not be used directly
# @test typeof(FAT.Simulations.load_snapshot(sim, sim.t[1], Val{:U})) == VectorField{2, Float64, typeof(mesh(sim))}
# @test typeof(FAT.Simulations.load_snapshot(sim, sim.t[1], Val{:p})) == ScalarField{2, Float64, typeof(mesh(sim))}

# # check loading snapshots with the sim[t] syntax. This should not be
# # used in performance critical code, because it is not type stable, as
# # the output depends on the value of the second argument
# @test typeof(sim[sim.t[1], :U]) == VectorField{2, Float64, typeof(mesh(sim))}
# @test typeof(sim[sim.t[1], :p]) == ScalarField{2, Float64, typeof(mesh(sim))}

# # check using integer and floats to load snapshots
# @test sim[sim.t[1], :p].internalField == sim[1, :p].internalField
# @test sim[sim.t[1], :p].boundaryField == sim[1, :p].boundaryField

# # check absence of snapshot
# @test_throws SnapshotError sim[1.123, :p]
# @test_throws SnapshotError sim[0, :p]
# @test_throws SnapshotError sim[length(sim)+1, :p]

# # test constructors of FieldIterator objects. Note 
# # that we need to use Val{:fieldname} for performance reasons
# fp = fields(sim, Val{:p})
# fU = fields(sim, Val{:U})
# fp = fields(sim, Val{:p}, 2:10)
# fU = fields(sim, Val{:U}, 2:10)
# @test_throws SnapshotError fields(sim, Val{:U}, 1:2:11)

# # test FieldIterator object
# @test eltype(fields(sim, Val{:U})) <: VectorField
# @test eltype(fields(sim, Val{:p})) <: ScalarField
# @test size(fields(sim, Val{:U}, 1:2)) == length(1:2)
# @test size(fields(sim, Val{:U}, 1:2:9)) == length(1:2:9)
# @test length(fields(sim, Val{:U}, 1:2)) == length(1:2)
# @test length(fields(sim, Val{:U}, 1:2:9)) == length(1:2:9)
# @test length(collect(fields(sim, Val{:U}, 1:2))) == length(1:2)
# @test length(collect(fields(sim, Val{:U}, 1:2:9))) == length(1:2:9)
# @test fields(sim, Val{:p}, 1:2:9)[2].internalField == sim[3.0, :p].internalField

# # test inference
# foo(sim, symb) = eltype(collect(fields(sim, Val{symb})))
# @test foo(sim, :U) <: VectorField
# @test foo(sim, :p) <: ScalarField
