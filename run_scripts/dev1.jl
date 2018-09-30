using CSV
using DataFrames
using Plots
using SeedMix
using Echidna

tablefile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Plant_data_updated_8_3_18.csv"
phylofile = "/Users/hawt0010/Projects/seed-mix/new_data_august_2018/Distance_ALL_7_31_18.csv"

sd = SpeciesData(tablefile, phylofile)
# objectivefuns = [get_phylo_dist, get_shannon, get_cost, get_consval, get_bloom]
# maxobjectives = [true, true, false, true, true]

objectivefuns = [get_phylo_dist, get_cost]
maxobjectives = [true, false]

nspecs = 15
mixweight = 10.0
evalfunction(g) = SeedMix.evaluate(objectivefuns, g, sd, nspecs, mixweight)

problem = Problem(length(objectivefuns), maxobjectives, 2*length(sd),
                  [MOGA_Real(0.0, 1.0) for i in 1:(2*length(sd))], evalfunction)
algo = NSGAII(problem, evalfunction, 100, 1000)

seedpop = init_pop(algo)
init_pd = [p.objectives[1] for p in seedpop]
init_cost = [p.objectives[2] for p in seedpop]

results = garun(algo, seedpop=seedpop)

post_pd = [p.objectives[1] for p in results]
post_cost = [p.objectives[2] for p in results]

gr()
scatter(init_cost, init_pd)
savefig("before_opt.png")

scatter(post_cost, post_pd)
savefig("after_opt.png")


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

all_selected = pool_spec_selections(results, nspecs, mixweight)
