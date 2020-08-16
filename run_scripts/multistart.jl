import Pkg
Pkg.activate("..")

using ArgParse
using Random
using SeedMixFunctions


function main()
    parsed_args = my_parse_cli_args()
    runid = randstring(6)
    println("Run string: $runid")

    project_root = ".."
    output_root = joinpath(project_root, "results/run_$(runid)")
    if ~isdir(output_root)
        mkpath(output_root)
    end
    logfile = joinpath(output_root, "log_$(runid).txt")

    mixreqs = MixRequirements(
        parsed_args["n"],
        parsed_args["W"],
        parsed_args["w"],
        parsed_args["c"])
    runparams = RunParams(
        parsed_args["t"],
        parsed_args["p"],
        parsed_args["g"])
    print_runspec(mixreqs, runparams)
    log_runspec(logfile, mixreqs, runparams)

    tablefile = joinpath(project_root, "data/species_data.csv")
    phylofile = joinpath(project_root, "data/phylo_dist.csv")
    sd = SpeciesData(tablefile, phylofile)

    multi_start(sd, mixreqs, runparams, output_root)

end


"Writes mixreqs and runparams to filename"
function log_runspec(filepath::String, mixreqs::MixRequirements, runparams::RunParams)
    open(filepath, "w") do f
        write(f, "MixRequirements\n")
        write(f, "    nspecs: $(mixreqs.nspecs)\n")
        write(f, "    totalweight: $(mixreqs.totalweight)\n")
        write(f, "    minweight: $(mixreqs.minweight)\n")
        write(f, "    maxcost: $(mixreqs.maxcost)\n")
        write(f, "RunParams\n")
        write(f, "    nstarts: $(runparams.nstarts)\n")
        write(f, "    popsize: $(runparams.popsize)\n")
        write(f, "    ngens: $(runparams.ngens)\n")
    end
end


"Prints mixreqs and runparams to screen"
function print_runspec(mixreqs::MixRequirements, runparams::RunParams)
    println("MixRequirements")
    println("    nspecs: $(mixreqs.nspecs)")
    println("    totalweight: $(mixreqs.totalweight)")
    println("    minweight: $(mixreqs.minweight)")
    println("    maxcost: $(mixreqs.maxcost)")
    println("RunParams")
    println("    nstarts: $(runparams.nstarts)")
    println("    popsize: $(runparams.popsize)")
    println("    ngens: $(runparams.ngens)")
end


"Parses commandline args"
function my_parse_cli_args()
    argsettings = ArgParseSettings()

    @add_arg_table argsettings begin
        "-n"
            help = "number of species in the mix"
            arg_type = Int64
            default = 25
        "-W"
            help = "total weight of the mix"
            arg_type = Float64
            default = 10.0
        "-w"
            help = "minimum weight of each species"
            arg_type = Float64
            default = 1.0/16.0
        "-t"
            help = "number of times to run the optimization"
            arg_type = Int64
            default = 1
        "-p"
            help = "number of indivs in the EA"
            arg_type = Int64
            default = 100
        "-g"
            help = "number of generations in the EA"
            arg_type = Int64
            default = 100
        "-c"
            help = "cost constraint (0.0 -> no constraint)"
            arg_type = Float64
            default = 0.0
    end

    return ArgParse.parse_args(ARGS, argsettings)
end


main()
