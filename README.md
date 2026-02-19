# Bruhat.jl

Tools to compute data related to reductive groups over finite fields. We heavily rely on Jean Michel's port of the computer algebra system [Chevie](https://github.com/jmichel7/Chevie.jl). 

Tools offered:

- `intersections_rational(G)` 

Returns the orders of Bruhat cells intersected with rational unipotent conjugacy classes. That is, returns the orders of $`B(\mathbb{F}_q) w B(\mathbb{F}_q) \cap C_g`$ where $`[w]\subseteq W`$ is a conjugacy class (with representative of minimal length among its conjugacy class) and $`C_g`$ is the $`G(\mathbb{F}_q)`$-conjugacy class of a unipotent element $`g\in G(\mathbb{F}_q)`$. This calculation follows [[Geck 2011](https://www.jstor.org/stable/43997674), Section 3].

- `intersections_geometric(G)` 

As above, but summing the columns according to which unipotent elements $`g\in G(\mathbb{F}_q)`$ lie in the same geometric conjugacy class. 


## Getting started

Download and install [Julia](https://julialang.org/downloads/). In the REPL (Julia's interactive command-line), copy-paste and run the below:

```julia
using Pkg; Pkg.add(url="https://github.com/baileywhitbread/Bruhat.jl")
```

This will install the `Bruhat.jl` package and its dependencies. To load the package, copy-paste and run the below:

```julia
using Bruhat
```

## Examples

### Intersection numbers when $G=G_2$

```julia
julia> using Bruhat

julia> G=coxgroup(:G,2)
G₂

julia> intersections_rational(G)

Rows labelled by conj classes C∈ConjCl(W)
Columns labelled by unipotent classes C⊆G(Fq)

┌───────┬─────────────────────────────────────────────────────────────────┐
│|BwB∩C|│1   A₁         Ã₁           G₂(a₁)    G₂(a₁)₍₂₁₎ G₂(a₁)₍₃₎     G₂│
├───────┼─────────────────────────────────────────────────────────────────┤
│A₀     │1 Φ₁Φ₃ (2q+1)q²Φ₁    (4q+1)q²Φ₁²/6 (2q+1)q²Φ₁²/2 q²Φ₁²Φ₂/3  q⁴Φ₁²│
│Ã₁     │0    0       q⁴Φ₁          q⁴Φ₁²/2       q⁴Φ₁²/2         0  q⁵Φ₁²│
│A₁     │0 q³Φ₁      q³Φ₁²     (q-2)q³Φ₁²/6       q⁴Φ₁²/2 q³Φ₁²Φ₂/3  q⁵Φ₁²│
│G₂     │0    0          0                0             0         0  q⁶Φ₁²│
│A₂     │0    0          0          q⁶Φ₁²/3             0  2q⁶Φ₁²/3  q⁸Φ₁²│
│A₁+Ã₁  │0    0      q⁶Φ₁² (q²-2q-2)q⁶Φ₁²/6  (q+2)q⁷Φ₁²/2   q⁶Φ₁⁴/3 q¹⁰Φ₁²│
└───────┴─────────────────────────────────────────────────────────────────┘

julia> intersections_geometric(G)

Rows labelled by conj classes C∈ConjCl(W)
Columns labelled by unipotent classes C⊆G(Fq)

┌───────┬────────────────────────────────────┐
│|BwB∩C|│1   A₁         Ã₁      G₂(a₁)     G₂│
├───────┼────────────────────────────────────┤
│A₀     │1 Φ₁Φ₃ (2q+1)q²Φ₁ (2q+1)q²Φ₁²  q⁴Φ₁²│
│Ã₁     │0    0       q⁴Φ₁       q⁴Φ₁²  q⁵Φ₁²│
│A₁     │0 q³Φ₁      q³Φ₁²       q⁴Φ₁²  q⁵Φ₁²│
│G₂     │0    0          0           0  q⁶Φ₁²│
│A₂     │0    0          0       q⁶Φ₁²  q⁸Φ₁²│
│A₁+Ã₁  │0    0      q⁶Φ₁²       q⁸Φ₁² q¹⁰Φ₁²│
└───────┴────────────────────────────────────┘
```

