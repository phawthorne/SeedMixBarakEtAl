using Documenter, SeedMix
push!(LOAD_PATH, "../src/")
makedocs(
    format = Documenter.HTML(),
    modules = [SeedMix],
    sitename = "SeedMix.jl",
    pages = ["Home" => "index.md"]
)
