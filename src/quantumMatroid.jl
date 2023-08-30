using Oscar
using Polymake
using Combinatorics
using AbstractAlgebra
using AbstractAlgebra.Generic

include("./matroid_relations.jl")
include("./utils.jl")
include("./multiSetMatroid.jl")

function partition_by(func::Function, sets)
    result = Dict()
    for set in sets
        key = func(set)
        if haskey(result,key)
            push!(result[key],set)
        else
            result[key] = [set]
        end
    end
    return result
end


function getRelations_rank(M::MultiSetMatroid)
    matroid = M.classic
    n = length(matroid)

    b =  [[] for _ in 1:n]
    nb =  [[] for _ in 1:n]
    rels= [[] for _ in 1:n]

    sts = powerset(matroid_groundset(matroid))   
    sts_by_rank = partition_by(x->rank(M,x),sts)
    sizes = unique(map(x -> length(x), sts))
    for size in 1:n
        tempGrdSet = reduce(vcat,[grdSet for i in 1:n])
        powerSet = unique(sort.(powerset(tempGrdSet,size,size)))

        prt = partition_by(x->rank(u,x),powerSet)
        for rnk in keys(prt)
            has_rank = prt[rnk]
            has_not_rank = setdiff(sts,has_rank)
        end

        for set in has_rank
            append!(b[size],collect(permutations(set)))
        end
        for nonset in has_not_rank
            append!(nb[size],collect(permutations(nonset)))
        end

    end
end




function getRelations(M::MultiSetMatroid,structure::Symbol=:bases)
    matroid = M.classic
    n = length(matroid)
    grdSet = matroid_groundset(matroid)

    b =  [[] for _ in 1:n]
    nb =  [[] for _ in 1:n]
    rels= [[] for _ in 1:n]

    sets  = eval(structure)(M)
    sizes = unique(map(x -> length(x), sets))
    for size in sizes 
        tempGrdSet = reduce(vcat,[grdSet for i in 1:n])
        powerSet = unique(sort.(powerset(tempGrdSet,size,size)))

        setsOfSize = filter(x->length(x)==size,sets)
        nonSets = setdiff(powerSet,setsOfSize) 

        for set in setsOfSize
            append!(b[size],collect(permutations(set)))
        end
        for nonset in nonSets
            append!(nb[size],collect(permutations(nonset)))
        end
        for set in b[size]
            for nonset in nb[size]
                rel = []
                for i in 1:size
                    push!(rel,(set[i],nonset[i]))
                end
                push!(rels[size],rel)
                rel = []
                for i in 1:size
                    push!(rel,(nonset[i],set[i]))
                end
                push!(rels[size],rel)
            end
        end

    end
    rels = unique.(rels)
    return reduce(vcat,rels)
end


function getIdeal(M::MultiSetMatroid,structure::Symbol=:bases; reduce::Bool=false)
    n = length(M.classic)  
    relations = getRelations(M,structure)
    result = matroid_relations(relations,n)
    reduce && interreduce!(result[3])
    Ideal = Oscar.FreeAssAlgIdeal(result[1],result[3])
    return result[1], Ideal, result[3], result[4]
end

getIdeal(M::Matroid,structure::Symbol=:bases; reduce::Bool=false)= getIdeal(MultiSetMatroid(M),structure,reduce=reduce)

function mygens(a)
  return [a.gens[Val(:O), i] for i in 1:ngens(a)]
end



function isInIdeal(ele::FreeAssAlgElem{T},gens::Vector{FreeAssAlgElem{T}},aut::AhoCorasickAutomaton) where T
   nrm = normal_form(ele,gens,aut) 
   return iszero(nrm) 
end

function isInIdeal(gens1::Vector{FreeAssAlgElem{T}},gens2::Vector{FreeAssAlgElem{T}},aut2::AhoCorasickAutomaton) where T
    for ele in gens1 
        isInIdeal(ele,gens2,aut2) || return false
    end
    return true
end
function generateSameIdeal(gens1::Vector{FreeAssAlgElem{T}},gens2::Vector{FreeAssAlgElem{T}}, aut1::AhoCorasickAutomaton, aut2::AhoCorasickAutomaton) where T
    isInIdeal(gens1,gens2,aut2) || return false
    isInIdeal(gens2,gens1,aut1) || return false
    return true
end

    
function isCommutative(gens,aut)
    Alg = parent(gens[1])
    n = ngens(Alg)
    for i in 1:n-1, j in i+1:n
        #print statement every 10%

        tst = Alg[i]*Alg[j] - Alg[j]*Alg[i]
        if !isInIdeal(tst,gens,aut)
            println("The Ideal is not commutative") 
            return tst 
        end
    end
    return true
end

function pwrSet(grdSet::Vector{<:Integer})
    powerSet_complete = Vector{Vector{Integer}}()
    for i in 1:length(grdSet)
        append!(powerSet_complete,multiset_combinations(grdSet,fill(i,length(grdSet)),i))
    end
    return powerSet_complete
end


function AdjacencyMatrix(v1::Vector{Vector{Integer}},v2::Vector{Vector{Integer}},m::Matroid)
    #this can be sped up with actual datastructures
    A = fill(1,length(v1),length(v2))
    circs =  circuits(m)
    for i in 1:length(v1), j in 1:length(v2)
        if !isempty(intersect(Set(v1[i]),Set(v2[j])))
            A[i,j] = 0
        else
            unio=union(v1[i],v2[j])
            for x in circs
                if issubset(x,unio) 
                    A[i,j] = 0
                    break
                end
            end
            
        end
    end        
    return A
end


#=
fan = fano_matroid()
nfan = non_fano_matroid()

Alg,Ideal_F,gens_F,aut_F = getIdeal(fan,:bases);
_, Ideal_NF, gens_NF, aut_NF = getIdeal(nfan,:bases);

res = isCommutative(Alg,gens_F,aut_F)
res = isCommutative(Alg,gens_NF,aut_NF)


isInIdeal(gens_NF,gens_F,aut_F)
isInIdeal(gens_F,gens_NF,aut_NF)

=#





