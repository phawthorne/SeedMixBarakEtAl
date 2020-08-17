import Pkg
Pkg.activate("..")

using Documenter, SeedMixFunctions
push!(LOAD_PATH, "../src/")
makedocs(
    format = Documenter.HTML(),
    modules = [SeedMixFunctions],
    sitename = "Documentation - Barak et al",
    pages = ["Home" => "index.md"]
)
