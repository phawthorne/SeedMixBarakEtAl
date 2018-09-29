using CSV
using DataFrames
using Params

@with_kw struct SpeciesData
    df::DataFrame
    phylo::Array{Float64, 2}
end

function SpeciesData(tablefile::String, phylofile::String)
    SpeciesData(
        CSV.File(tablefile) |> DataFrame
        zeros(5, 5)
    )
end


@with_kw struct MixData
    nspecs::Int64
    weight::Float64
    specindices::Vector{Int64}
    specweights::Vector{Float64}
    objectives::Vector{Float64}
end
