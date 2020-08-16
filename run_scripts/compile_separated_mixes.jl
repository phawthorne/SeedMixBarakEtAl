using CSV
using DataFrames
using DelimitedFiles
using SeedMixFunctions

"I should probably take a folder name as an argument"
function main()
    project_root = "/Users/hawt0010/Projects/julia-dev/SeedMix"
    tablefile = joinpath(project_root, "data/species_data.csv")
    phylofile = joinpath(project_root, "data/phylo_dist.csv")
    sd = SpeciesData(tablefile, phylofile)

    run_dir = "/Users/hawt0010/Projects/julia-dev/SeedMix/results/remote/maxcost-2500-2018-10-19"
    for (wroot, wdirs, wfiles) in walkdir(run_dir)
        if "mixes" in wdirs
            compile_mix_folder_nocsv(joinpath(wroot, "mixes"), sd)
        end
    end
end

"CSV.write is incredibly slow!"
function compile_mix_folder(foldername, sd)
    @show foldername
    files = readdir(foldername)
    specnames = sd.specpool
    df = DataFrame(species=specnames)

    for f in files
        @show f
        if f[1:3] !== "mix"
            continue
        end
        mixnum = parse(Int64, f[4:(end-4)])
        mixdf = CSV.read(joinpath(foldername, f))
        weight_vec = zeros(length(specnames))
        for (i, w) in zip(mixdf[:specindices], mixdf[:weights])
            weight_vec[i] = w
        end
        df[Symbol("mix$mixnum")] = weight_vec
    end

    CSV.write(joinpath(foldername, "compiled_mixes.csv"), df)
end

"This version is very fast!"
function compile_mix_folder_nocsv(foldername, sd)
    @show foldername
    files = readdir(foldername)
    mixfiles = [f for f in files if f[1:3]=="mix"]
    nmixes = length(mixfiles)

    specnames = sd.specpool
    nspecs = length(specnames)

    wmat = zeros(nspecs, nmixes)

    for i in 1:nmixes
        mixfile = joinpath(foldername, "mix$i.csv")
        mixdf = CSV.read(joinpath(foldername, mixfile))
        for (s, w) in zip(mixdf[:specindices], mixdf[:weights])
            wmat[s, i] = w
        end
    end

    open(joinpath(foldername, "compiled_mixes.csv"), "w") do io
        # header
        write(io, "species,")
        join(io, ("mix$i" for i in 1:nmixes), ",")
        write(io, "\n")
        #data
        for s in 1:nspecs
            write(io, "$(specnames[s]),") # first col
            join(io, wmat[s,:], ",")      # data rows
            write(io, "\n")
        end
    end
end

main()
