using Base.Test

using FAT.Simulation
using FAT.Meshes

for path in ("ldc_test_ascii", "ldc_test_binary")
    # use a test run
    sim = SimulationData(path, (1, 2))

    # check times
    @test times(sim) == collect(0:10)

    # check we return a mesh
    @test typeof(mesh(sim)) == Mesh{Float64}

    # check path
    @test casedir(sim) == path
end