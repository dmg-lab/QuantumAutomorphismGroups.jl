Pkg.activate(".")

include("./quantumMatroid.jl") # maybe use inclide("./src/quantumMatroid.jl") instead


M = uniform_matroid(3, 4)


Alg, Ideal, gens, aut = getIdeal(M);


isInIdeal(gens,gens,aut)

isCommutative(parent(gens[1]),gens,aut)

