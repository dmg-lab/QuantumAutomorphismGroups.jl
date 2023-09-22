using Oscar
using Combinatorics
using Oscar.AbstractAlgebra
using Oscar.AbstractAlgebra.Generic: FreeAssAlgElem

include("./utils.jl")
include("./multiSetMatroid.jl")
include("./save.jl")













function isInIdeal(ele::FreeAssAlgElem{T}, I::Oscar.FreeAssAlgIdeal{FreeAssAlgElem{T}},n::Int=3) where T<:FieldElem
        return ideal_membership(ele,I,n)
end

function isInIdeal(ele::FreeAssAlgElem{T}, gb::Vector{FreeAssAlgElem{T}}) where T<:FieldElem
    nrm = normal_form(ele,gb)
    return iszero(nrm) 
end

isInIdeal(ele :: FreeAssAlgElem{T}, gb::Vector{AbstractAlgebra.FreeAssAlgElem}) where T<:FieldElem = isInIdeal(ele, Vector{FreeAssAlgElem{QQFieldElem}}(gb))


function isInIdeal(gens1::Vector{FreeAssAlgElem{T}},gens2::Vector{FreeAssAlgElem{T}}, alt::Bool=false, n::Int=3) where T <:FieldElem
    if alt
        I=AbstractAlgebra.ideal(parent(gens2[1]),gens2)
        foreach(ele->isInIdeal(ele,I,n) || return false, gens1)        
    else
        for ele in gens1 
            isInIdeal(ele,gens2) || return false
        end
    end
    return true
end
function generateSameIdeal(gens1::Vector{FreeAssAlgElem{T}}, gens2::Vector{FreeAssAlgElem{T}}; alt::Bool=false, n::Int=3) where T <:FieldElem
    isInIdeal(gens1, gens2, alt, n) || return false
    isInIdeal(gens2, gens1, alt, n) || return false
    return true
end






function isCommutative(I::Union{Oscar.FreeAssAlgIdeal{FreeAssAlgElem{T}},Vector{FreeAssAlgElem{T}}}, all::Bool = true) where T <: FieldElem
    if typeof(I) == Oscar.FreeAssAlgIdeal{FreeAssAlgElem{T}}
        Alg = parent(gens(I)[1])
    else
        Alg = parent(I[1])
    end
    m = ngens(Alg)
    size = Int(sqrt(m))

    nCDict = Dict{FreeAssAlgElem,Bool}() # true means is in ideal, false means not in ideal
    c = true

    for i in 1:m-1, j in i+1:m
        tst = Alg[i]*Alg[j] - Alg[j]*Alg[i]

        j % size == i % size && continue
        div(i-1,size) == div(j-1,size) && continue

        if !isInIdeal(tst,I)
            nCDict[tst] = false
            all || return false, nCDict
            c = false
            continue
        end
        nCDict[tst] = true
    end
    return c, nCDict
end


function add_new_relation(relations::Vector{FreeAssAlgElem{T}}, new_relation::FreeAssAlgElem{T}) where T
    relations = copy(relations)
    push!(relations,new_relation)
    gb = AbstractAlgebra.Generic.groebner_basis(relations)
    println("Added new relation")
    return gb
end


function addMatroidRelations( 
    relationsToAdd::Vector{FreeAssAlgElem{T}} = FreeAssAlgElem{T}[],
    relations::Vector{FreeAssAlgElem{T}} = FreeAssAlgElem{T}[];
    path::String = "./MatroidRelationComputation") where T


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
            println("There is only $(freemem)% of free memory left. Saving... ")
            save(path * ".gb", (relations,addLater))
            if freemem <10
                println("Not enough memory to restart the computation. Please restart manually")
                exit()
            end
        end 
        
        relations = addMatroidRelations(addLater,relations)
    end
    #delete tmpfile
    if isfile(path * ".gb")
        rm(path * ".gb")
    end
    return relations
end




#=
gns , u , A = getMatroidRelations(fano_matroid())
gns
AbstractAlgebra.groebner_basis(gns)
gb = addMatroidRelations(gns)


restartComputation("./MatroidRelationComputation")


gns
isInIdeal(gns,gns)
=#


function restart()
    startup = """
        Base.ACTIVE_PROJECT[]=$(repr(Base.ACTIVE_PROJECT[]))
        Base.HOME_PROJECT[]=$(repr(Base.HOME_PROJECT[]))
        cd($(repr(pwd())))
        """
    cmd = `$(Base.julia_cmd()) -ie $startup`
    atexit(()->run(cmd;t=false))
    exit(0)
end






function restart(path::String)
    startup = """
        Base.ACTIVE_PROJECT[]=$(repr(Base.ACTIVE_PROJECT[]))
        Base.HOME_PROJECT[]=$(repr(Base.HOME_PROJECT[]))
        cd($(repr(pwd())))
        include($(repr("./quantumMatroid.jl")))
        println("Restarting Computation of $(path)")
        relations, relationsToAdd = load("$(path)" * ".gb")
        addMatroidRelations(relationsToAdd,relations)
        """
    cmd = `$(Base.julia_cmd()) -ie $startup`
    atexit(()->run(cmd))
    exit(0)
end

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

getMatroidRelations(M::Matroid,structure::Symbol=:bases, interreduce::Bool=false)= getMatroidRelations(MultiSetMatroid(M),structure,interreduce)

function isCommutative(M::Matroid,structure::Symbol=:bases, alt::Bool=false, n::Int=3)
    relsToAdd, _, A = getMatroidRelations(M,structure)
    if alt
        I = Oscar.ideal(A,relsToAdd)
        return isCommutative(I)
    end
    gb = AbstractAlgebra.groebner_basis(relsToAdd)
    return isCommutative(gb)
end


function getInfo(M::Matroid,alt::Bool=false, n::Int=3,save::Bool=false,path::String="./info.json")
    name =   String(M.pm_matroid.REVLEX_BASIS_ENCODING)
    aut_b,_ = isCommutative(M,:bases,alt,n)
    aut_c,_ = isCommutative(M,:circuits,alt,n)
    class_aut = automorphism_group(M)
    isTrans =  is_transitive(class_aut) 
    isSimple = is_simple(M)
    if save
        save(path,(name, Bool(aut_b), Bool(aut_c), isSimple, Oscar.describe(class_aut), isTrans))
    end

    return name, aut_b, aut_c, isSimple, Oscar.describe(class_aut), isTrans
end


#=
db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]

cursor=Polymake.Polydb.find(collection, Dict("RANK" => 2,"N_ELEMENTS"=>4))
Droids=Matroid.(cursor)








getInfo(Droids[1],true)
isCommutative(Droids[1],:flats,true)
String(Droids[1].pm_matroid.REVLEX_BASIS_ENCODING)

map(x->isCommutative(x,:rank,true)[1],Droids)

relsToAdd, gens, automat, u, A = getMatroidRelations(Droids[23])

isCommutativeSecondVersion(gens,3) #true 

=#



#gns , _ , _ = getMatroidRelations(fano_matroid())
#t = typeof(gns)
#
#gb = AbstractAlgebra.groebner_basis(gns)
##gb = addMatroidRelations(gns,t(),"./fano_backup")
#save("./fano_computation_alt.gb",gb)




