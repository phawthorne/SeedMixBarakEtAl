# Code documentation

SeedMixFunctions.jl is the included package that provides the mix evaluation
and GA run initiation functions. Use of these functions isn't required to 
run an analysis - the file `multistart.jl` in the scripts folder (documented [here](usage.html))
is the main way to run the analysis - but the functions are documented here
for comprehensiveness. 

## Data loading function
```@docs
SpeciesData(tablefile::String, phlofile::String)
```

## Data types
These are the structs that you need to run these analyses.
```@docs
SpeciesData
MixData
MixRequirements
```

## GA functions
Functions to run or post-process runs of the GA
```@docs
RunParams
multi_start
do_run
save_results
pool_spec_selections
```

## Evaluation functions
These are functions that evaluate a mix for a specific metric.

```@docs
get_phylo_dist(mix::MixData, sd::SpeciesData)
get_shannon
get_cost
get_consval
get_bloom
get_grass_spec_frac
get_grass_weight_frac
```

## Index

```@index
```
