using Oscar


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
    
Oscar.flats(M::MultiSetMatroid) = flats(M.classic)



#=
M = MultiSetMatroid(Oscar.fano_matroid())
rank(M.classic)
flats(M.classic)
flats(M)
closure(M,[1,2,3,3,3,3])

=#
