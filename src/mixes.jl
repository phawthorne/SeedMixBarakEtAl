using CSV
using DataFrames
using Parameters

@with_kw struct SpeciesData
    specpool::Vector{String}
    cost::Vector{Float64}
    conservatism::Vector{Float64}
    grass::Vector{Float64}
    bloom::Array{Float64, 2}
    phylo::Array{Float64, 2}
end

function SpeciesData(tablefile::String, phylofile::String)
    table = CSV.read(tablefile)
    phylo = CSV.read(phylofile)

    tablespecs = Set(table[:AcceptedName])
    phylospecs = Set(phylo[:AcceptedName])

    specpool = sort(table[:AcceptedName])

    tablelookup = Dict( s => findfirst(p->p==s, table[:AcceptedName]) for s in specpool )
    phylolookup = Dict( s => findfirst(p->p==s, phylo[:AcceptedName]) for s in specpool )

    cost = [table[tablelookup[s], :OneLbCost] for s in specpool]
    conservatism = [table[tablelookup[s], :Conservatism] for s in specpool]
    grass = [table[tablelookup[s], :Grass_Forb]=="G" ? 1.0 : 0.0 for s in specpool]

    months = [:Apr, :May, :Jun, :Jul, :Aug, :Sep, :Oct, :Nov]
    bloommat = zeros(length(specpool), length(months))
    for (i, s) in enumerate(specpool)
        si = tablelookup[s]
        for (j, m) in enumerate(months)
            bloommat[i, j] = ismissing(table[si, m]) ? 0.0 : 1.0
        end
    end

    phylomat = zeros(length(specpool), length(specpool))
    for (i, s) in enumerate(specpool)
        si = phylolookup[s]
        for (j, t) in enumerate(specpool)
            phylomat[i,j] = phylo[si, Symbol(t)] # have to do it this way because we know
                                             # the row order, but not col order
        end
    end

    SpeciesData(
        specpool,
        cost,
        conservatism,
        grass,
        bloommat,
        phylomat
    )
end


@with_kw struct MixData
    nspecs::Int64
    totalweight::Float64
    specindices::Vector{Int64}
    specweights::Vector{Float64}
    objectives::Vector{Float64}
end


function evaluate!(mix::MixData, sd::SpeciesData)
    pd = get_phylo_dist(mix, sd)
    shannon = get_shannon(mix, sd)
    cost = get_cost(mix, sd)
    consval = get_consval(mix, sd)
    bloom = get_bloom(mix, sd)
    grass_spec_frac = get_grass_spec_frac(mix, sd)
    grass_weight_frac = get_grass_weight_frac(mix, sd)
end
