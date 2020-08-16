# Seed Mix Optimization
Code repository for the optimization method described in Barak et al (2021). A 
[demonstration visualization and tool](https://phawthorne.github.io/computational-seed-mix-design/)
is also available.

## Installation
The core algorithms are written in the [Julia language](https://julialang.org/). Steps to running the code:

1. Download Julia (code was tested with Julia 1.5).
2. Clone this code repository
3. Set up the Julia environment
    - Open Julia
    - Enter these commands:
    ```
    julia> include Pkg
    julia> activate .
    julia> instantiate
    ```
    - Close Julia
4. Set up Julia to be able to call it from the command line. 
5. Navigate to the scripts directory.
6. Run command with desired arguments (see [documentation]() for options)
    ```
    % julia multistart.jl [args]
    ```

## Other notes
This code depends on [Echidna](https://github.com/phawthorne/Echidna), which implements the underlying 
NSGA-ii optimization algorithm, and will be installed automatically when the environment is instantiated.