using CSV
using DataFrames
using DelimitedFiles
using Echidna


struct RunParams
    nstarts::Int64
    popsize::Int64
    ngens::Int64
end


"Perform a multi-start run of the GA"
function multi_start(sd::SpeciesData, mixreqs::MixRequirements,
    runparams::RunParams, output_root::String)
    #= Do the runs =#
    results = Vector{Vector{Solution}}()
    for i in 1:runparams.nstarts
        println("GA start: $i")
        runresult = do_run(sd, mixreqs, runparams)
        save_results(joinpath(output_root, "run$i"), runresult, mixreqs, sd)
        push!(results, runresult)
    end

    #= Write results =#
    pooled_results = Archive(compare_pareto_dominance, Vector{Solution}())
    for (i,r) in enumerate(results)
        insert_solutions!(pooled_results, r)
        @show length(pooled_results.solutions)
    end
    save_results(joinpath(output_root, "pooled"),
                 pooled_results.solutions, mixreqs, sd)
end


"Perform a single run of the GA with given data and parameters"
function do_run(sd::SpeciesData, mixreqs::MixRequirements, runparams::RunParams)
    @unpack popsize, ngens = runparams

    objectivefuns = [get_cost, get_phylo_dist, get_bloom, get_shannon, get_consval]
    maxobjectives = [false, true, true, true, true]
    has_cost_constraint = mixreqs.maxcost > 0.0

    evalfunction(g) = SeedMixFunctions.evaluate(objectivefuns, g, sd, mixreqs)

    problem = Problem(
        length(objectivefuns),
        maxobjectives,
        2*length(sd),
        [MOGA_Real(0.0, 1.0) for i in 1:(2*length(sd))],
        evalfunction,
        has_cost_constraint)
    hall_of_fame = Archive(compare_pareto_dominance, Vector{Solution}())
    archive_frequency = 100000
    algo = NSGAII(problem, evalfunction, popsize, ngens, hall_of_fame, archive_frequency)
    seedpop = init_pop(algo)
    results = garun(algo, seedpop=seedpop)

    return results
end


"Performs union of species found in results' mixes"
function pool_spec_selections(results, mixreqs::MixRequirements)
    specs = Set{Int64}()
    for r in results
        m = genome_to_mix(r.x, mixreqs)
        for s in m.specindices
            push!(specs, s)
        end
    end
    return specs
end


"Writes output files for genomes, objectives, and mix details"
function save_results(output_folder::String, results::Vector{Solution},
                      mixreqs::MixRequirements, specdata::SpeciesData)
    if !isdir(output_folder)
        mkpath(output_folder)
    end

    # save matrix of genome vals
    genomes = zeros(length(results), length(results[1].x))
    for (i, r) in enumerate(results)
        genomes[i,:] = r.x
    end
    writedlm(joinpath(output_folder, "genomes.csv"), genomes, ",")

    mixes = [genome_to_mix(r.x, mixreqs) for r in results]

    # save table of objective values
    indiv = collect(1:length(results))
    generation = [r.generation for r in results]
    cost = [r.objectives[1] for r in results]
    phylo_dist = [r.objectives[2] for r in results]
    bloom = [r.objectives[3] for r in results]
    shannon = [r.objectives[4] for r in results]
    consval = [r.objectives[5] for r in results]
    grass_spec_frac = [get_grass_spec_frac(m, specdata) for m in mixes]
    grass_weight_frac = [get_grass_weight_frac(m, specdata) for m in mixes]
    df = DataFrame(indiv=indiv, generation=generation,
                   cost=cost, phylo_dist=phylo_dist, bloom=bloom,
                   shannon=shannon, consval=consval,
                   grass_spec_frac=grass_spec_frac,
                   grass_weight_frac=grass_weight_frac)
    CSV.write(joinpath(output_folder, "objectives.csv"), df)

    # save per-mix seed and weight tables
    if !isdir(joinpath(output_folder, "mixes"))
        mkpath(joinpath(output_folder, "mixes"))
    end
    specs_per_mix = mixreqs.nspecs
    weight_per_mix = mixreqs.totalweight
    speclist = results[1].problem.eval_fn.sd.specpool
    for (i, r) in enumerate(results)
        filename = joinpath(output_folder, "mixes", "mix$i.csv")
        m = genome_to_mix(r.x, mixreqs)
        df = DataFrame(
            specindices = m.specindices,
            specnames = speclist[m.specindices],
            weights = m.specweights
        )
        CSV.write(filename, df)
    end
end
