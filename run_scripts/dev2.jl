using BenchmarkTools
using CSV
using DataFrames
using DelimitedFiles
using Parameters
using Random
# using Plots
using SeedMix
using Echidna

struct RunParams
    nstarts::Int64
    popsize::Int64
    niters::Int64
end

function main()
    runid = randstring(6)
    println("Run string: $runid")

    if length(ARGS) >= 1
        arg1 = parse(Int64, ARGS[1])
        println("Using input var: $arg1")
    else
        arg1 = 25
    end

    if "HOSTNAME" in keys(Base.EnvDict())
        project_root = "/home/hawt0010/hawt0010/Projects/SeedMix"
    else
        project_root = "/Users/hawt0010/Projects/julia-dev/SeedMix"
    end
    output_root = joinpath(project_root, "results/run_$(runid)")
    if ~ isdir(output_root)
        mkpath(output_root)
    end
    logfile = joinpath(output_root, "log_$(runid).txt")

    mixreqs = MixRequirements(arg1, 10.0, 1.0/16.0)
    runparams = RunParams(10, 100, 100)
    logrun(logfile, mixreqs, runparams)

    tablefile = joinpath(project_root, "data/species_data.csv")
    phylofile = joinpath(project_root, "data/phylo_dist.csv")
    sd = SpeciesData(tablefile, phylofile)

    #= Do the runs =#
    results = Vector{Vector{Solution}}()
    for i in 1:runparams.nstarts
        println("GA start: $i")
        runresult = dorun(sd, mixreqs, runparams)
        push!(results, runresult)
    end

    #= Write results =#
    pooled_results = Archive(compare_pareto_dominance, Vector{Solution}())
    for (i,r) in enumerate(results)
        save_results(joinpath(output_root, "run$i"), r, mixreqs)
        insert_solutions!(pooled_results, r)
        @show length(pooled_results.solutions)
    end
    save_results(joinpath(output_root, "pooled"), pooled_results.solutions, mixreqs)
end

function dorun(sd::SpeciesData, mixreqs::MixRequirements, runparams::RunParams)
    @unpack popsize, niters = runparams

    objectivefuns = [get_cost, get_phylo_dist, get_bloom]
    maxobjectives = [false, true, true]

    evalfunction(g) = SeedMix.evaluate(objectivefuns, g, sd, mixreqs)

    problem = Problem(length(objectivefuns), maxobjectives, 2*length(sd),
                      [MOGA_Real(0.0, 1.0) for i in 1:(2*length(sd))], evalfunction)
    hall_of_fame = Archive(compare_pareto_dominance, Vector{Solution}())
    archive_frequency = 100000
    algo = NSGAII(problem, evalfunction, popsize, niters, hall_of_fame, archive_frequency)
    seedpop = init_pop(algo)
    results = garun(algo, seedpop=seedpop)

    return results
end

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

function save_results(output_folder::String, results::Vector{Solution}, mixreqs::MixRequirements)
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
    cost = [r.objectives[1] for r in results]
    phylo_dist = [r.objectives[2] for r in results]
    bloom = [r.objectives[3] for r in results]
    df = DataFrame(indiv=indiv, cost=cost, phylo_dist=phylo_dist, bloom=bloom)
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

function logrun(filepath::String, mixreqs::MixRequirements, runparams::RunParams)
    open(filepath, "w") do f
        write(f, "MixRequirements\n")
        write(f, "\tnspecs: $(mixreqs.nspecs)\n")
        write(f, "\ttotalweight: $(mixreqs.totalweight)\n")
        write(f, "\tminweight: $(mixreqs.minweight)\n")
        write(f, "RunParams\n")
        write(f, "\tnstarts: $(runparams.nstarts)\n")
        write(f, "\tpopsize: $(runparams.popsize)\n")
        write(f, "\tniters: $(runparams.niters)\n")
    end
end



main()
