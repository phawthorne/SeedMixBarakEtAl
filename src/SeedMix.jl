module SeedMix

include("mixes.jl")
export SpeciesData, MixData, length,
    get_phylo_dist, get_shannon, get_cost, get_consval, get_bloom,
    genome_to_mix, evaluate, random_mix

end # module
