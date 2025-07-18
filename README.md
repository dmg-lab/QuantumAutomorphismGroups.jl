# QuantumAutomorphismGroups.jl

QuantumAutomorphismGroups.jl is part of an ongoing research project on quantum automorphisms of matroids [[1]](#1). It contains both the code needed to compute the quantum automorphism groups of matroids and the data we computed for several matroids. All the computations were done using the OSCAR system, and the data was encoded using OSCAR's internal machinery [[2]](#2).


## Usage

When you first start up, you need to activate the project folder and instantiate it.

```bash
julia --project /path/to/folder -e 'using Pkg; Pkg.instantiate()'
```

Most of its original functionality is now available in the OSCAR system. See the OSCAR documentation for more information: [Link](https://docs.oscar-system.org/stable/Combinatorics/matroids/#Quantum-Automorphisms).
However, it still serves as a database of computable quantum automorphism groups of matroids, accessible via the `load_dt` function.

### Accessing the Database

The following example shows how to load a list of all computed quantum automorphism groups:

```julia-repl
julia> using QuantumAutomorphismGroups, Oscar, DataFrames

julia> dt = load_dt();

julia> df2 = select(dt.data, ["Name", "Aut_bases_-1", "Aut_circuits_-1","Aut_circuits_6", "Aut_bases_6"]);

julia> df2 = filter(row -> row["Aut_bases_6"] !== missing || !ismissing(row["Aut_bases_-1 "]), df2)
250×5 DataFrame
 Row │ Name            Aut_bases_-1  Aut_circuits_-1  Aut_circuits_6  Aut_bases_6 
     │ String15        Bool?         Bool?            Bool?           Bool?       
─────┼────────────────────────────────────────────────────────────────────────────
   1 │ r1n7_01              missing          missing         missing        false
   2 │ r0n1_1                  true          missing            true      missing 
   3 │ r1n2_3                  true             true            true      missing 
   4 │ r1n3_1                  true             true            true      missing 
   5 │ r1n3_3                  true             true            true      missing 
   6 │ r1n3_7                  true          missing            true      missing 
   7 │ r1n4_01                 true             true            true      missing 
   8 │ r1n4_03                 true             true            true      missing 
   9 │ r1n4_07                 true             true            true      missing 
  10 │ r1n4_0f                false            false           false      missing 
  11 │ r1n5_01                false            false           false      missing 
  12 │ r1n5_03                 true             true            true      missing 
  13 │ r1n5_07                 true             true            true      missing 
  14 │ r1n5_0f                false            false           false      missing 
  15 │ r1n5_1f                false            false           false      missing 
  16 │ r1n6_01                false            false           false      missing 
  17 │ r1n6_03                false            false           false      missing 
  18 │ r1n6_07                 true             true            true      missing 
  19 │ r1n6_0f                false            false           false      missing 
  ⋮  │       ⋮              ⋮               ⋮               ⋮              ⋮
 233 │ r3n7_07ffffffe       missing          missing            true        false
 234 │ r3n7_07fffffff         false          missing            true        false
 235 │ r3n7_0965b96ef       missing          missing         missing         true
 236 │ r3n7_1bcffbffb       missing          missing            true         true
 237 │ r3n7_1beffbfff       missing          missing            true         true
 238 │ r3n7_3f7eefd6f       missing          missing         missing         true
 239 │ r3n7_3f7eefd7f       missing          missing            true         true
 240 │ r3n7_3f7eefdff       missing          missing            true         true
 241 │ r3n7_3f7eeffff       missing          missing         missing         true
 242 │ r3n7_3f7effdfe       missing          missing            true         true
 243 │ r3n7_3f7effdff       missing          missing            true         true
 244 │ r3n7_3f7efffbf       missing          missing            true         true
 245 │ r3n7_3f7efffff       missing          missing            true         true
 246 │ r3n7_3f7fffff7       missing          missing            true         true
 247 │ r3n7_3f7ffffff       missing          missing            true         true
 248 │ r3n7_3ffff7fff       missing          missing            true         true
 249 │ r3n7_3ffffffff       missing          missing            true         true
 250 │ r3n7_7ffffffff       missing          missing         missing         true
```

To access the matroid by its name, you can use the OSCAR function `matroid_from_matroid_hex':

```julia-repl
julia> matroid_from_matroid_hex("r1n4_01")
Matroid of rank 1 on 4 elements
```

To return the complete data on a matroid use the function `load_dict`:

```julia-repl
julia> y = load_dict(uniform_matroid(1,3))
Dict{String, Any} with 4 entries:
  "Aut_bases_-1_timed"   => 0.123875
  "Aut_circuits_6_timed" => 1.34087
  "Aut_bases_-1"         => FreeAssociativeAlgebraElem[u[1, 1]^2 - u[1, 1], u[1, 1]*u[1,…
  "Aut_circuits_6"       => FreeAssociativeAlgebraElem[u[1, 3] + u[2, 3] + u[3, 3] - 1, …
```
## References

<a id="1">[1]</a>
Daniel Corey, Michael Joswig, Julien Schanz, Marcel Wack, Moritz Weber </br>
Quantum automorphisms of matroids </br>
[arXiv](https://arxiv.org/abs/2312.13464)

<a id="2">[2]</a>
Antony Della Vecchia, Michael Joswig and Benjamin Lorenz </br>
A FAIR File Format for Mathematical Software (2023) </br>
[arXiv](https://arxiv.org/abs/2309.00465)

