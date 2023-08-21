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

   #= 
u = uniform_matroid(3,4)
sts = powerset(matroid_groundset(u))
collect(sts)
prt = partition_by(x->rank(u,x),sts)
#merge all sets that dont have a specific rank
rnk = 3
has_rank = prt[rnk]
has_not_rank = setdiff(sts,has_rank)
for key in keys(prt)
    prt[key] = reduce(vcat,prt[key])
end
prt
=#

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

#=
M = MultiSetMatroid(fano_matroid())
rels = getRelations(M,:bases)
=#

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

function generateSameIdeal(I1,I2, n)
    gb1 = mygens(I1)
    gb2 = mygens(I2)
    for ele in gb1 
        if !ideal_membership(ele,I2,n)
            return false
        end
    end
    for ele in gb2 
        if !ideal_membership(ele,I1,n)
            return false
        end
    end
    return true
end

function isInIdeal(I1,I2, n)
    gb1 = mygens(I1)
    gb2 = mygens(I2)
    for ele in gb1 
        if !ideal_membership(ele,I2,n)
            return false
        end
    end
    return true
end



function isCommutative(Alg,I,n)
    for i in 1:ngens(Alg),j in 1:ngens(Alg)
        tst = Alg[i]*Alg[j] - Alg[j]*Alg[i]
        ideal_membership(tst,I,n) || return false
    end
    return true
end

#= Uniform matroid rank 2 with 3 elements

uni = uniform_matroid(4,5)

Alg, Ideal, gens, aut = getIdeal(uni,:bases);
Alg, Ideal_C,gb_C, AHO_C = getIdeal(uni,:circuits);
x= gb_C[end-20]

res = normal_form_with_rep(x,gb_B,AHO_B)
normal_form_with_rep(res[2],gb_B,AHO_B)

ideal_membership(res[2],Ideal_B,3)
=#


matroid = uniform_matroid(3,7)
grdSet = matroid.groundset
c = 5
d = 3

powerSet_c = multiset_combinations(grdSet,fill(c,length(grdSet)),c)
s_c = length(powerSet_c)





powerSet_d = multiset_combinations(grdSet,fill(d,length(grdSet)),d)
s_d = length(powerSet_d)


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



function getRelations2(M::MultiSetMatroid,structure::Symbol=:bases)
    matroid = M.classic
    grdSet = matroid.groundset
end

#=
uni = uniform_matroid(4,5)
multiUni = MultiSetMatroid(uni) 



Alg, Ideal_B,gb_B, AHO_B = getIdeal(uni,:bases);
tst = 1

for x in Alg
    tst*=x
end

Alg, Ideal_C,gb_C, AHO_C = getIdeal(uni,:circuits);
gb_B
x = gb_C[end-20]
x= Alg[1]*Alg[7]*Alg[10]*Alg[13]
x1 = Alg[1]*Alg[8]-Alg[8]*Alg[1]

GB =  Generic.groebner_basis(gb_B)
auts = Generic.AhoCorasickAutomaton(Vector{Int}[])

rels = []
rel_count = 0
for i in 1:length(GB)
    rel_count += 1
    add_new_relation!!(GB,auts,GB[i],rel_count)
end

res = normal_form_with_rep(x,GB,auts)
res[1]
normal_form_with_rep(GB[end],gb_B,AHO_B)

x1 
res[2]
sol = normal_form_with_rep(GB[4142],vcat(gb_B,GB),AHO_B)
sol[2]
sol[1]


uni = uniform_matroid(3,4)
multiUni = MultiSetMatroid(uni) 


FreeAssAlg, Ideal_B,gb_B, AHO_B = getIdeal(uni,:bases);
typeof(Ideal_B)
rat = FreeAssAlg.base_ring

typeof(FreeAssAlg.base_ring)
typeof(FreeAssAlg.base_ring)

FreeAssAlg
save("test",rat)
typeof( Ideal_B.gens[1])


typeof(Ideal_B)
save("tst.tst",gb_B[1])


Alg, Ideal_C,gb_C, AHO_C = getIdeal(uni,:circuits);
x= Alg[1]*Alg[7]*Alg[10]*Alg[13]
su = gb_C[116]
x1 = x * (gb_C[116]+1)

arr = collect(AbstractAlgebra.monomials(x1))
for y in arr
    println(ideal_membership(y,Ideal_B,5))
end
import AbstractAlgebra: QQ, elem_type

tst = Alg[1]*Alg[7]*Alg[10]*Alg[13]*Alg[10]


GB =  Generic.groebner_basis(gb_B)
auts = Generic.AhoCorasickAutomaton(Vector{Int}[])

rels = []
rel_count = 0
for i in 1:length(GB)
    rel_count += 1
    add_new_relation!!(GB,auts,GB[i],rel_count)
end
GB
res = normal_form_with_rep(x,GB,auts)






GB = ideal_membership(tst2,Ideal_B,4)

println.(gb_C[10:120])
GB =  Generic.groebner_basis(gb_B)

sol = normal_form_with_rep(x1,GB,AHO_B)
generateSameIdeal(Ideal_C,Ideal_B,4)

idel_membership(x1,Ideal_B,5)
=#







#=
gens = ["x","y"]
A, g = FreeAssociativeAlgebra(Oscar.QQ,gens)
f = g[1]^2+g[2]
g = g[2]

x = Oscar.ideal(A,[f,g])
typeof(x) 
R, x = QQ["x"]
typeof(x^2)
=#





