using Oscar: interreduce
using Oscar
using Polymake
using Combinatorics
using AbstractAlgebra
using AbstractAlgebra.Generic

include("./matroid_relations.jl")
include("./utils.jl")
include("./multiSetMatroid.jl")
include("./save.jl")




getIdeal(M::Matroid,structure::Symbol=:bases; reduce::Bool=false)= getIdeal(MultiSetMatroid(M),structure,reduce=reduce)



function isInIdeal(ele::FreeAssAlgElem{T}, gb::Vector{FreeAssAlgElem{T}}, alt::Bool=false, n::Int=3) where T<:FieldElem
    if alt
        return ideal_membership(ele,gb,n)
    end
    nrm = normal_form(ele,gb)
    return iszero(nrm) 
end

function isInIdeal(gens1::Vector{FreeAssAlgElem{T}},gens2::Vector{FreeAssAlgElem{T}}, alt::Bool=false, n::Int=3) where T <:FieldElem
    for ele in gens1 
        isInIdeal(ele,gens2,alt,n) || return false
    end
    return true
end
function generateSameIdeal(gens1::Vector{FreeAssAlgElem{T}}, gens2::Vector{FreeAssAlgElem{T}}; alt::Bool=false, n::Int=3) where T <:FieldElem
    isInIdeal(gens1, gens2, alt, n) || return false
    isInIdeal(gens2, gens1, alt, n) || return false
    return true
end



function isCommutative(gb::Vector{FreeAssAlgElem{T}}, alt::Bool=false, n::Int=3) where T <: FieldElem

    Alg = parent(gb[1])
    m = ngens(Alg)
    
    for i in 1:m-1, j in i+1:m
        
        tst = Alg[i]*Alg[j] - Alg[j]*Alg[i]
        if !isInIdeal(tst,gb,alt,n)
            println("The Ideal is not commutative") 
            return tst 
        end
    end
    return true
end






function add_new_relation(relations::Vector{FreeAssAlgElem{T}}, new_relation::FreeAssAlgElem{T}) where T
    relations = copy(relations)
    push!(relations,new_relation)
    return AbstractAlgebra.Generic.groebner_basis(relations)
end


function addMatroidRelations( 
    relationsToAdd::Vector{FreeAssAlgElem{T}} = FreeAssAlgElem{T}[],
    relations::Vector{FreeAssAlgElem{T}} = FreeAssAlgElem{T}[]) where T

    tmpfileString = "./MatroidRelationComputation"

    #Get the first Quarter of relations_indices
    quarter = div(length(relationsToAdd),4)
    if quarter == 0
        quarter = length(relationsToAdd)
    end

    addNow = relationsToAdd[1:quarter]
    addLater = relationsToAdd[quarter+1:end]

    for rel in addNow
        relations = add_new_relation(relations, rel)
    end


    if addLater == []
        return relations
    else
        freemem = round(getFreePercentageOfMemory() * 100,digits=2)
        if freemem < 40 
            println("There is only $(freemem)% of free memory left. Saving to $(tmpfileString * ".gb")")

            save(tmpfileString * ".gb", (relations,addLater))
            println("Saved with $(string(length(addNow))) new relations added.\n There are now $(string(length(addLater))) relations left")
            if freemem <10
                println("Not enough memory to restart the computation. Please restart manually")
                exit()
            end

        end 
        
        relations = addMatroidRelations(addLater,relations)
    end

    return relations
end


#=
gns, automat = getQuantumPermutationGroup(3)
gns
interreduce!(gns)

gb = addMatroidRelations(gns)
restartComputation("./MatroidRelationComputation")


gns
isInIdeal(gns,gns)
=#





function restartComputation(path::String)
    println("Restarting Computation of $(path)")
    relations, relationsToAdd = load(path * ".gb")
    addMatroidRelations(relationsToAdd,relations)
    return

end
    





function getMatroidRelations(
    M::MultiSetMatroid,
    structure::Symbol=:bases,
    interreduce::Bool=false)
    
    relation_indices = getRelations(M,structure)
    relation_transformed,  u, A = getQuantumPermutationGroup(length(M.classic),interreduce)

    for relation in relation_indices
        temp = one(A)
        for gen in relation
            temp = temp * u[gen[1], gen[2]]
        end
        push!(relation_transformed,temp)
    end
  
    return relation_transformed,  u, A

end

getMatroidRelations(M::Matroid,structure::Symbol=:bases)= getMatroidRelations(MultiSetMatroid(M),structure)

#=
relsToAdd, u, A = getMatroidRelations(uniform_matroid(2,3))
relsToAdd_C, u_C, A_C = getMatroidRelations(uniform_matroid(2,3),:circuits)

gb = AbstractAlgebra.groebner_basis(relsToAdd)
gb_C = AbstractAlgebra.groebner_basis(relsToAdd_C)

isInIdeal(gb,gb)
isCommutative(gb) # true

isInIdeal(gb_C,gb_C)
isCommutative(gb_C) # true

generateSameIdeal(gb,gb_C) # true
=#

#=
relsToAdd, u, A = getMatroidRelations(uniform_matroid(2,4))
relsToAdd_C, u_C, A_C = getMatroidRelations(uniform_matroid(2,4),:circuits)
relsToAdd_F, u_F, A_F = getMatroidRelations(uniform_matroid(2,4),:flats)

gb = AbstractAlgebra.groebner_basis(relsToAdd)
gb_C = AbstractAlgebra.groebner_basis(relsToAdd_C)
gb_F = AbstractAlgebra.groebner_basis(relsToAdd_F)


isInIdeal(gb,gb)
isCommutative(gb) # true

isInIdeal(gb_C,gb_C)
isCommutative(gb_C) # true

generateSameIdeal(gb,gb_C) # true
=#
#=

relsToAdd, u, A = getMatroidRelations(uniform_matroid(3,4))
relsToAdd_C, u_C, A_C = getMatroidRelations(uniform_matroid(3,4),:circuits)
relsToAdd_F, u_F, A_F = getMatroidRelations(uniform_matroid(3,4),:flats)

isCommutativeSecondVersion(relsToAdd,3) #true
isCommutativeSecondVersion(relsToAdd_C,3) #true
isCommutativeSecondVersion(relsToAdd_F,3) #true
=#
function isCommutative(M::Matroid,structure::Symbol=:bases, alt::Bool=false, n::Int=3)
    relsToAdd, u, A = getMatroidRelations(M,structure)
    gb = AbstractAlgebra.groebner_basis(relsToAdd)
    return isCommutative(gb)
end


function getNameAndCommutative(M::Matroid)
    name =   String(M.pm_matroid.REVLEX_BASIS_ENCODING)
    aut_b = isCommutative(M,:bases)
    aut_c = isCommutative(M,:circuits)
    return name, aut_b, aut_c
end




#=
db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]

cursor=Polymake.Polydb.find(collection, Dict("RANK" => 3,"SIMPLE"=>false,"N_ELEMENTS"=>7))
Droids=Matroid.(cursor)

getNameAndCommutative_alt(Droids[1])
String(Droids[1].pm_matroid.REVLEX_BASIS_ENCODING)
isCommutative_alt(Droids[1],:bases)

map(getNameAndCommutative_alt,Droids)

relsToAdd, gens, automat, u, A = getMatroidRelations(Droids[23])

isCommutativeSecondVersion(gens,3) #true 

=#


gns, automat = getMatroidRelations(MultiSetMatroid(fano_matroid()),:bases)
gb = addMatroidRelations(gns)

