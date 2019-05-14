# RepresentativeDaysFinders.jl


A Julia Module to  ...

## Getting started

### Installation

From the julia REPL, run

```julia
v1.1> add https://git.vito.be/scm/emark-epdst/representativedaysfinder.jl.git
```

This will install `representativedaysfinder.jl` as well as all its dependencies

### Upgrading

To upgrade to the most recent version of `representativedaysfinder.jl`, run


```julia
v1.1> up
```

Alternatively, you can specify a branch, as in

```julia
v1.1> up
```


### Usage

Run:

```julia
using RepresentativeDaysFinders
using JuMP
using GLPK
##################################################################################
# Specify location of config-file
##################################################################################
config_file = normpath(joinpath(@__DIR__, "scenarios", "DE_DK_2015_1.yaml"))

findRepresentativeDays(config_file, with_optimizer(GLPK.Optimizer; presolve=true, msg_lev=GLPK.MSG_ALL))
```
## Trouble shooting
If issues with GR-engine occur just build GR package:

```julia

build GR

```


## Documentation

...
## Reporting Issues and Contributing

...

## License

...
