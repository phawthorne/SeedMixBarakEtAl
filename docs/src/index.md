# Seed Mix Code

Documentation for `SeedMix.jl`, a julia package for optimizing and analyzing seed mixes for prairie restoration.

```@contents
```

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
