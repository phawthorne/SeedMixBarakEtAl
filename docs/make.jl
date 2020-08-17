import Pkg
Pkg.activate("..")

using Documenter, SeedMixFunctions
push!(LOAD_PATH, "../src/")
makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [SeedMixFunctions],
    sitename = "Documentation - Barak et al",
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "Code documentation" => "SeedMixFunctions.md"
    ]
)
