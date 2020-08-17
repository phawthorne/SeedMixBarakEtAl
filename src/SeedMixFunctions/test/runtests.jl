# import Pkg
# Pkg.activate("..")

using Test, SeedMixFunctions

println("Testing")

mix = MixData(10, 10.0,
    collect(range(1, length=10)),
    ones(10))

@test length(mix) == 10
