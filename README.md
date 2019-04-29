# RepresentativeDaysFinders.jl


A Julia Module to  

## Getting started

### Installation

While in the Julia REPL, press the "]" key to get to the package manager REPL mode (this is similar to `using Pkg`).

```julia
(v1.1) pkg> add git@gitlab.mech.kuleuven.be:UCM/representativedaysfinder.jl.git
```

This will install `representativedaysfinder.jl` as well as all its dependencies

For devoloping the package clone it to <path_to_git_clone> and

```julia
(v1.1) pkg> dev <path_to_git_clone>
```

### Upgrading

To upgrade to the most recent version of `representativedaysfinder.jl`, run


```julia
(v1.1) pkg> up RepresentativeDaysFinders
```

### Usage

Run:

```julia
using RepresentativeDaysFinders
using JuMP
using GLPK

# Specify location of config-file

config_file = normpath(joinpath(@__DIR__, "scenarios", "DE_DK_2015_1.yaml"))

findRepresentativeDays(config_file, with_optimizer(GLPK.Optimizer; presolve=true, msg_lev=GLPK.MSG_ALL))
```
## Trouble shooting
If issues with GR-engine occur just build GR package:

```julia
(v1.1) pkg> build GR
```

## Developers

```julia
julia> cd("<path_to_git_clone>")
(v1.1) pkg> activate <path_to_git_clone>
```
If you add packages now they're added to the `Project.toml` file.

## Documentation

...

## Reporting Issues and Contributing

If you have any issues, report them on GitLab and mention '@steffen.kaminski' (will speed up response to the issue).

## License

...
