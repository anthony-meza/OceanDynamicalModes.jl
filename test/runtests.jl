using OceanDynamicalModes
using Test

@testset "OceanDynamicalModes.jl" begin
    # Write your tests here.
    dz = 0.05
    depth = collect(0 + dz:dz:1. - 0.05)
    Nsq = ones(length(depth))
    nmodes = 3
    wmodes, pmodes, ce = dynmodes(Nsq, depth, dz, nmodes)
end
