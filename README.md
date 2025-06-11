# Bruhat.jl

Some tools to compute data related to reductive groups over finite fields. For instance:

- Sizes of Bruhat cells intersected with unipotent conjugacy classes: $|B(\mathbb{F}_q) w B(\mathbb{F}_q) \cap G(\mathbb{F}_q)\cdot g|$.


## Getting started

Download and install [Julia](https://julialang.org/downloads/). In the REPL (Julia's interactive command-line), copy-paste and run the below:

```julia
using Pkg; Pkg.add(url="https://github.com/baileywhitbread/Bruhat.jl")
```

This will install the `Bruhat.jl` package and its dependencies. To load the package, copy-paste and run the below:

```julia
using Bruhat
```