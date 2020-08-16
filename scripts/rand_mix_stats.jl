using CSV
using DataFrames
using Statistics
using Random
using SeedMixFunctions
using Echidna

tablefile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Plant_data_updated_8_3_18.csv"
phylofile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Distance_ALL_7_31_18.csv"

sd = SpeciesData(tablefile, phylofile)
objectivefuns = [get_phylo_dist, get_shannon, get_cost, get_consval, get_bloom]
maxobjectives = [true, true, false, true, true]
mixweight = 10.0

function make_pop(specs_in_mix, n_mixes)
    return [[Random.rand(MOGA_Real(0.0, 1.0)) for i in 1:(2*length(sd))] for i in 1:10]
end

function find_stats(min_specs, max_specs, n_mixes)
    means = zeros(max_specs, 5)
    stddevs = zeros(max_specs, 5)
    medians = zeros(max_specs, 5)
    for specs_in_mix in min_specs:max_specs
        sample_indivs = [[Random.rand(MOGA_Real(0.0, 1.0)) for i in 1:(2*length(sd))] for i in 1:n_mixes]
        mr = MixRequirements(specs_in_mix, mixweight, 1/16)
        evalfunction(g) = SeedMix.evaluate(objectivefuns, g, sd, mr)
        scores = [evalfunction(g) for g in sample_indivs]
        means[specs_in_mix, :] .= [mean([s[i] for s in scores]) for i in 1:5]
        stddevs[specs_in_mix, :] .= [std([s[i] for s in scores]) for i in 1:5]
        medians[specs_in_mix, :] .= [median([s[i] for s in scores]) for i in 1:5]
    end
    return (means, stddevs, medians)
end

function write_file(filepath, stats, min_specs, max_specs)
    open(filepath, "w") do f
        write(f, "SpecsInMix,PD,Shannon,Cost,ConsVal,Bloom\n")
        for sp in min_specs:max_specs
            write(f, "$(sp),$(stats[sp,1]),$(stats[sp,2]),$(stats[sp,3]),$(stats[sp,4]),$(stats[sp,5])\n")
        end
    end
end

function run_analysis(min_specs::Int, max_specs::Int, n_mixes::Int)
    means, stddevs, medians = find_stats(min_specs, max_specs, n_mixes)
    write_file("random_samples/means.csv", means, 10, 100)
    write_file("random_samples/stddevs.csv", stddevs, 10, 100)
    write_file("random_samples/medians.csv", medians, 10, 100)
end

min_specs = 10
max_specs = 100
n_mixes = 100000
run_analysis(min_specs, max_specs, n_mixes)
