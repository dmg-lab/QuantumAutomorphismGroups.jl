using Oscar
using Oscar.AbstractAlgebra
using Oscar.AbstractAlgebra.Generic

include("./utils.jl")
include("./multiSetMatroid.jl")
include("./save.jl")

export isInIdeal,
    generateSameIdeal,
    isCommutative,
    isCommutativeExtra,
    getMatroidRelations



function isInIdeal(ele::FreeAssAlgElem{T}, I::Oscar.FreeAssAlgIdeal{<:FreeAssAlgElem{T}},n::Int=3) where T<:FieldElem

        return ideal_membership(ele,I,n)
end

function isInIdeal(ele::FreeAssAlgElem{T}, gb::Vector{<:FreeAssAlgElem{T}}) where T<:FieldElem
    nrm = normal_form(ele,gb)
    return iszero(nrm) 
end



function isInIdeal(gens1::Vector{<:FreeAssAlgElem{T}}, gens2::Vector{<:FreeAssAlgElem{T}}, alt::Bool=false, n::Int=3) where T <:FieldElem
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

function generateSameIdeal(gens1::Vector{<:FreeAssAlgElem{T}}, gens2::Vector{<:FreeAssAlgElem{T}}; alt::Bool=false, n::Int=3) where T <:FieldElem
    isInIdeal(gens1, gens2, alt, n) || return false
    isInIdeal(gens2, gens1, alt, n) || return false
    return true
end


function isCommutative(I::Vector{<:FreeAssAlgElem{T}}, all::Bool = true) where T <: FieldElem

    Alg = parent(I[1])
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

function isCommutative(M::Matroid, structure::Symbol = :bases, alt::Bool = false, n::Int = 3)
    relsToAdd, _, A = getMatroidRelations(M,structure)
    if alt
        I = Oscar.ideal(A,relsToAdd)
        return isCommutative(I)
    end
    gb = AbstractAlgebra.groebner_basis(relsToAdd)
    return isCommutative(gb)
end

function isCommutativeExtra(I::Vector{<:FreeAssAlgElem{T}}, all::Bool = false) where T <: FieldElem
    Alg = parent(I[1])
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
            tst2 = Alg[i]*Alg[j] - Alg[j]*Alg[i]*Alg[j]
            if !isInIdeal(tst2,I)
                nCDict[tst2] = false
                all || return false, nCDict
                c = false
                continue
            end
        end
        nCDict[tst] = true
    end
    return c, nCDict
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

getMatroidRelations(M::Matroid, structure::Symbol=:bases, interreduce::Bool=false)= getMatroidRelations(MultiSetMatroid(M),structure,interreduce)

function getInfo(M::Matroid, alt::Bool=false, n::Int=3, save::Bool=false, path::String="./info.json")
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
gns , u , A = getMatroidRelations(uniform_matroid(2,4))
typeof(gns)
isInIdeal(gns,gns)
isCommutative(gns)
isCommutativeExtra(gns)
getInfo(uniform_matroid(1,2))


=#
