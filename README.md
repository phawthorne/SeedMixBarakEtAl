# Seed Mix Optimization

[![Build Status](https://travis-ci.com/phawthorne/SeedMixBarakEtAl.svg?branch=master)](https://travis-ci.com/phawthorne/SeedMixBarakEtAl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://phawthorne.github.io/SeedMixBarakEtAl/dev)

Code repository for the optimization method described in Barak et al (2021). A 
[demonstration visualization tool](https://github.com/phawthorne/computational-seed-mix-design/)
is also available.

## Installation
The core algorithms are written in the [Julia language](https://julialang.org/). Steps to running the code:

1. Download Julia (code was tested with Julia 1.5).
1. Set up Julia to be able to call it from the command line.
    - On Mac, for example, add `alias julia="/Applications/Julia-1.5.app/Contents/Resources/julia/bin/julia"` to your `.zshrc` file.
1. Clone this code repository using git.
1. Set up the Julia environment
    - Navigate to this repo's folder in the command line
    - Open Julia
    - Enter these commands:
    ```
    julia> include Pkg
    julia> Pkg.activate(".")
    julia> Pkg.instantiate()
    julia> exit()
    ```
1. Now navigate to the scripts directory.
1. Run the optimization script (`julia multistart.jl`) with desired arguments 
(see [documentation]() for options)

## Other notes
This code depends on [Echidna](https://github.com/phawthorne/Echidna), which implements the underlying 
NSGA-ii optimization algorithm, and will be installed automatically when the environment is instantiated.

Other files in the scripts folder, `rand_mix_stats.jl` and `rand_mix_stats_plots.py` were used to generate
the z-scores for calculating normalized/adjusted bloom coverage objective values. These will require some
updates to run on other computers. Also note that `rand_mix_stats_plots.py` requires Python 3
with Pandas, Numpy, and Matplotlib installed.

