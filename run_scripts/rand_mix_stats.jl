include("../src/SeedMix.jl")
using CSV
using DataFrames
using Statistics
using Random
using .SeedMix
using Echidna

tablefile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Plant_data_updated_8_3_18.csv"
phylofile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Distance_ALL_7_31_18.csv"

sd = SpeciesData(tablefile, phylofile)
objectivefuns = [get_phylo_dist, get_shannon, get_cost, get_consval, get_bloom]
maxobjectives = [true, true, false, true, true]
mixweight = 10.0

function find_stats(specs_in_mix::Int, n_mixes::Int)
    sample_indivs = [[Random.rand(MOGA_Real(0.0, 1.0)) for i in 1:(2*length(sd))] for i in 1:10]
    mr = MixRequirements(specs_in_mix, mixweight, 1/16)
    evalfunction(g) = SeedMix.evaluate(objectivefuns, g, sd, mr)
    scores =     [evalfunction(g) for g in sample_indivs]
    medians =    [median([s[i] for s in scores]) for i in 1:5]
    return medians
end

function run_analysis(filepath::String, min_specs::Int, max_specs::Int, n_mixes::Int)
    open(filepath, "w") do f
        write(f, "SpecsInMix,PD,Shannon,Cost,ConsVal,Bloom\n")
        for sp in min_specs:max_specs
            stats = find_stats(sp, n_mixes)
            write(f, "$(sp),$(stats[1]),$(stats[2]),$(stats[3]),$(stats[4]),$(stats[5])\n")
        end
    end
end

fp = "test.csv"
min_specs = 10
max_specs = 100
n_mixes = 100000
run_analysis(fp, min_specs, max_specs, n_mixes)
