# Usage

The main way to run the code is through the `scripts/multistart.jl` script, which 
loads the appropriate data and executes an optimization run. The script must
be called through the command line from the `scripts` folder:

```
% cd scripts
% julia multistart.jl
```

This will initiate a run of the optimization with default arguments. In order to customize
the analysis, additional command line arguments can be specified.

## multistart.jl arguments
The available arguments are:
```
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
```

The analyses in the paper were:
```
% julia multistart.jl -t 20 -c 2500 -n 25 -p 200 -g 5000
% julia multistart.jl -t 20 -c 2500 -n 32 -p 200 -g 5000
% julia multistart.jl -t 20 -c 2500 -n 40 -p 200 -g 5000
% julia multistart.jl -t 20 -c 5000 -n 25 -p 200 -g 5000
% julia multistart.jl -t 20 -c 5000 -n 32 -p 200 -g 5000
% julia multistart.jl -t 20 -c 5000 -n 40 -p 200 -g 5000
```

It does take a while to run with this many starts, generations, and individuals in the population,
so we recommend experimenting with the default before running this size of problem. 