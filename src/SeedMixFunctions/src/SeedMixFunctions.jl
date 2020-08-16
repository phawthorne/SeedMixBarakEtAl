module SeedMix

include("mixes.jl")
export SpeciesData, MixData, MixRequirements, length,
    get_phylo_dist, get_shannon, get_cost, get_consval, get_bloom,
    get_grass_spec_frac, get_grass_weight_frac,
    genome_to_mix, evaluate, random_mix

end # module
