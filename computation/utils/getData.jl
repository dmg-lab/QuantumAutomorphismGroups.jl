
using Oscar

include("../../src/quantumMatroid.jl")
include("../../src/multiSetMatroid.jl")





#=
M = uniform_matroid(3,4)
structure = :bases
=#

function computeQuantumAutoMorphism(M::MultiSetMatroid, structure::Symbol)
    duration = @elapsed result = begin 
        n = length(M.classic)  
        relations = getRelations(M,structure)
        matroid_relations(relations,n)
    end

    return result[1],result[3],result[4], duration
end

computeQuantumAutoMorphism(M::Matroid, structure::Symbol) = computeQuantumAutoMorphism(MultiSetMatroid(M), structure)

Alg , gens, aut, dur = computeQuantumAutoMorphism(M, structure)

