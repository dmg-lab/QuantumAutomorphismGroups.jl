using Oscar

export
    getQuantumPermutationGroup



@doc raw"""

    MultiSetMatroid(matroid::Matroid)

Construct a MultiSetMatroid from a Matroid.

# Examples
```julia
QuantumAutomorphismGroups.MultiSetMatroid(fano_matroid())

# output

QuantumAutomorphismGroups.MultiSetMatroid(Matroid of rank 3 on 7 elements)
```
"""
struct MultiSetMatroid
    classic::Oscar.Matroid
end

Oscar.length(M::MultiSetMatroid) = length(M.classic)
Oscar.bases(M::MultiSetMatroid) = bases(M.classic)
Oscar.independent_sets(M::MultiSetMatroid)= independent_sets(M.classic)

function Oscar.circuits(M::MultiSetMatroid)
    matroid = M.classic
    groundset = matroid_groundset(matroid)
    cc = Oscar.circuits(matroid)
    lops = Oscar.loops(matroid)
    for ele in 1:length(matroid)
        ele in lops && continue
        push!(cc,[ele,ele])
    end
    return cc
end

@doc raw"""
    rank(M::MultiSetMatroid)

Return the rank of the MultiSetMatroid `M`.

# Examples
```jldoctest
julia> rank(MultiSetMatroid(fano_matroid()))
3
```
"""
function Oscar.rank(M::MultiSetMatroid)
    return Oscar.rank(M.classic)
end
    

@doc raw"""
    rank(M::MultiSetMatroid,set::Vector{Int})

Return the rank of a Multiset `set` as part of the MultiSetMatroid `M`.

# Examples
```jldoctest
julia> rank(MultiSetMatroid(fano_matroid()),[1,1])
1
```
"""
function Oscar.rank(M::MultiSetMatroid,set::Vector{Int})
    return Oscar.rank(M.classic,unique(set))
end

function Oscar.closure(M::MultiSetMatroid,set::Vector{Int})
    return closure(M.classic,unique(set))
end

@doc raw"""
    flats(M::MultiSetMatroid)

Return the flats of the MultiSetMatroid `M`.

# Examples
```jldoctest 
M = MultiSetMatroid(Oscar.fano_matroid())
flats(M)

```
"""
Oscar.flats(M::MultiSetMatroid) = flats(M.classic)


@doc raw"""
    getRelations(M::MultiSetMatroid,structure::Symbol=:bases)

Get the indices of the relations that define the quantum automorphism group of a matroid for a given structure.

# Examples
```jldoctest
M = uniform_matroid(3,4)
idx = QuantumAutomorphismGroups.getRelations(M,:bases)
length(idx)

# output

1920
```
"""
function getRelations(M::MultiSetMatroid,structure::Symbol=:bases)
    structure == :rank && return getRelations_rank(M)

    matroid = M.classic
    n = length(matroid)
    grdSet = matroid_groundset(matroid)
    
    b =  [[] for _ in 1:n]
    nb =  [[] for _ in 1:n]
    rels= [[] for _ in 1:n]

    sets  = eval(structure)(M)
    sizes = unique(map(x -> length(x), sets))
    for size in sizes 
        size == 0 && continue
        tempGrdSet = reduce(vcat,[grdSet for i in 1:n])
        powerSet = unique(sort.(powerset(tempGrdSet,size,size)))

        setsOfSize = filter(x->length(x)==size,sets)
        nonSets = setdiff(powerSet,setsOfSize) 

        for set in setsOfSize
            append!(b[size],collect(Oscar.permutations(set)))
        end
        for nonset in nonSets
            append!(nb[size],collect(Oscar.permutations(nonset)))
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
    return Vector{Vector{Tuple{Int,Int}}}(reduce(vcat,rels))
end

getRelations(M::Matroid,structure::Symbol=:bases)=getRelations(MultiSetMatroid(M),structure)


function partition_by(func::Function, sets)
    result = Dict{Int,Vector{Vector{Int}}}()
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
    grdSet = matroid_groundset(matroid)
    b =  [[] for _ in 1:n]
    nb =  [[] for _ in 1:n]
    rels= [[] for _ in 1:n]

    sts = powerset(matroid_groundset(matroid))   
    sts_by_rank = partition_by(x->rank(M,Vector{Int}(x)),sts)
    for size in 1:n
        tempGrdSet = reduce(vcat,[grdSet for i in 1:n])
        powerSet = unique(sort.(powerset(tempGrdSet,size,size)))
        prt = partition_by(x->rank(M,x),powerSet)
        for rnk in keys(prt)
            has_rank = prt[rnk]
            has_not_rank = setdiff(powerSet,has_rank)

            for set in has_rank
                append!(b[size],collect(Oscar.permutations(set)))
            end
            for nonset in has_not_rank
                append!(nb[size],collect(Oscar.permutations(nonset)))
            end


        end
        for set in b[size]
            length(set) == 0 && continue
            for nonset in nb[size]
                length(nonset) == 0 && continue
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
    return Vector{Vector{Tuple{Int,Int}}}(reduce(vcat,rels))

end

getRelations_rank(M::Matroid)=getRelations_rank(MultiSetMatroid(M))


@doc raw"""
    getQuantumPermutationGroup(n::Int)

Get the relations that define the quantum permutation group on `n` elements.

# Examples
```jldoctest
rels,_, u, A = QuantumAutomorphismGroups.getQuantumPermutationGroup(3)
length(rels)
# output

20
```
"""
function getQuantumPermutationGroup(n::Int)
    generator_strings = String[]
    for i in 1:n, j in 1:n
            push!(generator_strings, "u[$i,$j]")
    end
    A, g = free_associative_algebra(Oscar.QQ, generator_strings)
    u = Matrix{elem_type(A)}(undef, n, n)
    for i in 1:n, j in 1:n
            u[i, j] = g[(i-1)*n + j]
    end
    rels_by_type = Dict{Symbol, Vector{elem_type(A)}}()
    rels_by_type[:idempotent] = elem_type(A)[]
    rels_by_type[:row_sum] = elem_type(A)[]
    rels_by_type[:col_sum] = elem_type(A)[]
    rels_by_type[:zero_divisor] = elem_type(A)[]
    #Idempotent relations
    for i in 1:n, j in 1:n
        new_relation = u[i, j] * u[i, j] - u[i, j]
        push!(rels_by_type[:idempotent], new_relation)
            for k in 1:n
                    if k != j
                        new_relation = u[i,j] * u[i, k]
                        push!(rels_by_type[:zero_divisor], new_relation)
                        new_relation = u[j, i]*u[k, i]
                        push!(rels_by_type[:zero_divisor], new_relation)
                    end
            end
    end

    #row and column sum relations
    for i in 1:n
        new_relation_row = -1
        new_relation_col = -1
        for k in 1:n
            new_relation_row += u[i,k]
            new_relation_col += u[k,i]
        end
        push!(rels_by_type[:row_sum], new_relation_row)
        push!(rels_by_type[:col_sum], new_relation_col)

    end
    relations = elem_type(A)[]
    for rel_type in keys(rels_by_type)
        relations = vcat(relations, rels_by_type[rel_type])
    end
    return Vector{elem_type(A)}(relations),rels_by_type, u, A
end
