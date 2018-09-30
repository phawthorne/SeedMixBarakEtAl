using CSV
using DataFrames
using Parameters
using Random
import Base.length

#= Species Data =#
@with_kw struct SpeciesData
    specpool::Vector{String}
    cost::Vector{Float64}
    conservatism::Vector{Float64}
    grass::Vector{Bool}
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
    grass = [table[tablelookup[s], :Grass_Forb]=="G" for s in specpool]

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

function length(sd::SpeciesData)
    return length(sd.specpool)
end

#= MixData =#
@with_kw struct MixData
    nspecs::Int64
    totalweight::Float64
    specindices::Vector{Int64}
    specweights::Vector{Float64}
end

function MixData(m::Vector{Tuple{Int64, Float64}})
    nspecs = length(m)
    si = [t[1] for t in m]
    sw = [t[2] for t in m]
    MixData(length(m), sum(sw), si, sw)
end

function length(md::MixData)
    return md.nspecs
end



#= Evaluation functions =#
function get_phylo_dist(mix::MixData, sd::SpeciesData)
    @unpack specindices = mix
    @unpack phylo = sd
    pd = 0.0
    for s in specindices
        for t in specindices
            pd += phylo[s,t]
        end
    end
    return pd /= 2.0
end

function get_shannon(mix::MixData, sd::SpeciesData)
    @unpack specweights = mix
    p = specweights ./ sum(specweights)
    return shannon = exp(-1.0 * sum(p .* log.(p)))
end

function get_cost(mix::MixData, sd::SpeciesData)
    @unpack specindices, specweights = mix
    @unpack cost = sd
    mixcost = 0.0
    for (s, w) in zip(specindices, specweights)
        mixcost += w * cost[s]
    end
    return mixcost
end

function get_consval(mix::MixData, sd::SpeciesData)
    @unpack specindices, specweights = mix
    @unpack conservatism = sd
    mixcons = 0.0
    for (s, w) in zip(specindices, specweights)
        mixcons =+ w * conservatism[s]
    end
    return mixcons
end

function get_bloom(mix::MixData, sd::SpeciesData, penaltyfun=sqrt)
    @unpack specindices, specweights = mix
    @unpack bloom = sd

    blooms_month_count = sum(bloom[specindices,:], dims=1)
    return sum(penaltyfun.(blooms_month_count))
end

function get_grass_spec_frac(mix::MixData, sd::SpeciesData)
    @unpack specindices = mix
    @unpack grass = sd

    mixG = Float64(sum(grass[specindices]))
    return mixG/mix.nspecs
end

function get_grass_weight_frac(mix::MixData, sd::SpeciesData)
    @unpack specindices, specweights = mix
    @unpack grass = sd

end

#= GA helper functions =#

function genome_to_mix(genome::Vector{Float64}, nspecs, mixweight)
    p = Int64(length(genome)/2)
    spec_priority = genome[1:p]
    spec_weight = genome[(p+1):end]

    spec_choice = sortperm(spec_priority, rev=true)[1:nspecs]
    sw = spec_weight[spec_choice]

    MixData(nspecs,
            mixweight,
            spec_choice,
            sw .* (mixweight/sum(sw)))
end

function evaluate(objectivefuns::Vector{Function}, mix::MixData, sd::SpeciesData)
    result = [f(mix, sd) for f in objectivefuns]
    return result
end

function evaluate(objectivefuns::Vector{Function}, genome::Vector{Float64},
                  sd::SpeciesData, nspecs::Int64, mixweight::Float64)
    mix = genome_to_mix(genome, nspecs, mixweight)
    result = [f(mix, sd) for f in objectivefuns]
    return result
end

function random_mix(sd::SpeciesData, nspecs::Int64, weight::Float64)
    specs = randperm(length(sd))[1:nspecs]
    w = rand(nspecs)
    w = w * (weight/sum(w))
    return MixData(nspecs, weight, specs, w)
end
