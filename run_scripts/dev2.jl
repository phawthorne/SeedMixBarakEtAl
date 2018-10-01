using BenchmarkTools
using CSV
using DataFrames
using DelimitedFiles
# using Plots
using SeedMix
using Echidna

function main()
    tablefile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Plant_data_updated_8_3_18.csv"
    phylofile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Distance_ALL_7_31_18.csv"

    sd = SpeciesData(tablefile, phylofile)

    nruns = 10
    nspecs = 25
    mixweight = 10.0
    popsize = 200
    niters = 1000

    pools = Vector{Set{Int64}}()
    results = Vector{Vector{Solution}}()
    for i in 1:nruns
        print("run: $i")
        runresult = dorun(sd, popsize, niters)
        push!(results, runresult)
        speclist = pool_spec_selections(runresult, nspecs, mixweight)
        push!(pools, speclist)
    end

    pooled_results = Archive(compare_pareto_dominance, Vector{Solution}())
    for (i,r) in enumerate(results)
        save_results("results/cost_vs_pd/run$i", r)
        insert_solutions!(pooled_results, r)
        @show length(pooled_results.solutions)
    end
    save_results("results/cost_vs_pd/pooled", pooled_results.solutions)

end

function dorun(sd::SpeciesData, popsize::Int64, niters::Int64)
    objectivefuns = [get_phylo_dist, get_cost]
    maxobjectives = [true, false]

    nspecs = 15
    mixweight = 10.0
    evalfunction(g) = SeedMix.evaluate(objectivefuns, g, sd, nspecs, mixweight)

    problem = Problem(length(objectivefuns), maxobjectives, 2*length(sd),
                      [MOGA_Real(0.0, 1.0) for i in 1:(2*length(sd))], evalfunction)
    hall_of_fame = Archive(compare_pareto_dominance, Vector{Solution}())
    archive_frequency = 100000
    algo = NSGAII(problem, evalfunction, popsize, niters, hall_of_fame, archive_frequency)
    seedpop = init_pop(algo)
    results = garun(algo, seedpop=seedpop)

    return results
end

function pool_spec_selections(results, nspecs, weight)
    specs = Set{Int64}()
    for r in results
        m = genome_to_mix(r.x, nspecs, weight)
        for s in m.specindices
            push!(specs, s)
        end
    end
    return specs
end

function save_results(output_folder::String, results::Vector{Solution})
    if !isdir(output_folder)
        mkpath(output_folder)
    end

    # save matrix of genome vals
    genomes = zeros(length(results), length(results[1].x))
    for (i, r) in enumerate(results)
        genomes[i,:] = r.x
    end
    writedlm(joinpath(output_folder, "genomes.csv"), genomes, ",")

    # save table of objective values
    indiv = collect(1:length(results))
    phylo_dist = [r.objectives[1] for r in results]
    cost = [r.objectives[2] for r in results]
    df = DataFrame(indiv=indiv, phylo_dist=phylo_dist, cost=cost)
    CSV.write(joinpath(output_folder, "objectives.csv"), df)

    # save per-mix seed and weight tables

    if !isdir(joinpath(output_folder, "mixes"))
        mkpath(joinpath(output_folder, "mixes"))
    end
    specs_per_mix = results[1].problem.eval_fn.nspecs
    weight_per_mix = results[1].problem.eval_fn.mixweight
    speclist = results[1].problem.eval_fn.sd.specpool
    for (i, r) in enumerate(results)
        filename = joinpath(output_folder, "mixes", "mix$i.csv")
        m = genome_to_mix(r.x, specs_per_mix, weight_per_mix)
        df = DataFrame(
            specindices = m.specindices,
            specnames = speclist[m.specindices],
            weights = m.specweights
        )
        CSV.write(filename, df)
    end

end


main()
