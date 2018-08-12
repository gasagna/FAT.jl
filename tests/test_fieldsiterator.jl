using Base.Test

using FAT.Simulation

for path in ("ldc_test_ascii", "ldc_test_binary")
    # simulation object
    sim = SimulationData(path, (1, 2))

    # construct iterator object for velocity skipping field at t=0
    us = vfields(sim, :U, 1:10)

    # check length
    @test length(us) == 10

    # check from data file ldc_test_ascii/1/U
    @test us[1][1].internalField[1] == 0.000259015397
    @test us[1][2].internalField[end] == -0.147795835
end