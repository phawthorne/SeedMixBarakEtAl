# Installation

The core algorithms are written in the [Julia language](https://julialang.org/). Steps to running the code:

+ Download Julia (code was tested with Julia 1.5).
+ Set up Julia to be able to call it from the command line.
    - On Mac, for example, add `alias julia="/Applications/Julia-1.5.app/Contents/Resources/julia/bin/julia"` to your `.zshrc` file.
+ Clone this code repository using git.
+ Set up the Julia environment
    - Navigate to this repo's folder in the terminal
    - Open Julia in the terminal
    - Enter these commands:
```
julia> include Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
julia> exit()
```
+ Now change to the scripts directory in the terminal.
+ Run the optimization script (`julia multistart.jl`) with desired arguments 
(see [Usage](@ref) for options)
