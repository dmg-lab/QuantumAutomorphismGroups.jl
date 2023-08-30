#!/usr/bin/julia
using Oscar
import AbstractAlgebra.Generic: FreeAssociativeAlgebra, FreeAssAlgElem, AhoCorasickAutomaton, insert_keyword!, normal_form


function add_new_relation!!(relations::Vector{AbstractAlgebra.Generic.FreeAssAlgElem{T}}, aut::AhoCorasickAutomaton, new_relation::FreeAssAlgElem{T}) where T
    normalized = normal_form(new_relation, relations, aut)
    if !iszero(normalized)
        push!(relations, normalized)
        insert_keyword!(aut, normalized.exps[1], length(relations))
    end
end
#Type richtig machen Vector{Tupe{Int,Int}}
function matroid_relations(relation_indices::Vector{Any}, n::Int)
    #Setup
    generator_strings = String[]
    relation_count = 0
    for i in 1:n, j in 1:n
            push!(generator_strings, "u[$i,$j]")
    end
    A, g = FreeAssociativeAlgebra(Oscar.QQ, generator_strings)
    u = Matrix{elem_type(A)}(undef, n, n)
    for i in 1:n, j in 1:n
            u[i, j] = g[(i-1)*n + j]
    end
    relations = elem_type(A)[]
    aut = AhoCorasickAutomaton(Vector{Int}[])

    #Squared relations
    for i in 1:n, j in 1:n
        new_relation = u[i, j] * u[i, j] - u[i, j]
        relation_count += 1
        if length(relations) == 0
            push!(relations, new_relation)
            insert_keyword!(aut, new_relation.exps[1], length(relations))
        else
            add_new_relation!!(relations, aut, new_relation, relation_count)
        end
            for k in 1:n
                    if k != j
                        relation_count += 1
                        new_relation = u[i,j] * u[i, k]
                        add_new_relation!!(relations, aut, new_relation, relation_count)
                        new_relation = u[j, i]*u[k, i]
                        relation_count += 1
                        add_new_relation!!(relations, aut, new_relation, relation_count)
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
        relation_count += 1
        add_new_relation!!(relations,aut,new_relation_row,relation_count)
        relation_count += 1
        add_new_relation!!(relations,aut,new_relation_col,relation_count)
    end

    #relations from the matroid
    m = length(relation_indices)
    relation_transformed = Vector{elem_type(A)}()
    for relation in relation_indices
        temp = one(A)
        for gen in relation
            temp = temp * u[gen[1], gen[2]]
        end
        push!(relation_transformed,temp)
    end
    
    #interreduce!(relation_transformed)

    i = 1
    for relation in relation_transformed
        i+=1
        if i % 1000 == 0
            GC.gc()
            percent = round(100*i/m,digits=2)
            println("i = $i, percent = $percent %")
        end
        relation_count += 1
        add_new_relation!!(relations, aut, relation, relation_count)
    end



    return A, u, relations, aut
    
end
function check_commutativity(u::Matrix{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}}, gb::Vector{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}}, aut::AhoCorasickAutomaton)
        for i in 1:size(u)[1]
                for j in 1:size(u)[2]
                        for k in 1:size(u)[1]
                                for l in 1:size(u)[2]
                                        if !iszero(normal_form(u[i, j]*u[k, l] - u[k, l] * u[i, j], gb, aut))
                                                return false
                                        end
                                end
                        end
                end
        end
        return true

end

function check_commutativity(u::Matrix{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}}, gb::Vector{AbstractAlgebra.Generic.FreeAssAlgElem{QQFieldElem}})
        for i in 1:size(u)[1]
                for j in 1:size(u)[2]
                        for k in 1:size(u)[1]
                                for l in 1:size(u)[2]
                                        if !iszero(normal_form(u[i, j]*u[k, l] - u[k, l] * u[i, j], gb))
                                                return false
                                        end
                                end
                        end
                end
        end
        return true
end


function normal_form_with_rep(
    f::FreeAssAlgElem{T},
    g::Vector{FreeAssAlgElem{T}},
    aut::AhoCorasickAutomaton,
) where {T}
rep_dict = Dict{Int, FreeAssAlgElem{T}}()
    R = parent(f)
    rexps = AbstractAlgebra.Generic.Monomial[]
    rcoeffs = T[]
    while length(f) > 0
        ok, left, right, match_index = AbstractAlgebra.Generic.gb_divides_leftmost(f.exps[1], aut)
        if ok
            qi = AbstractAlgebra.divexact_right(f.coeffs[1], g[match_index].coeffs[1])
            f = AbstractAlgebra.Generic._sub_rest(f, AbstractAlgebra.Generic.mul_term(qi, left, g[match_index], right), 1)
            rep_dict[match_index] = AbstractAlgebra.Generic.mul_term(qi, left, g[match_index], right)
        else
            push!(rcoeffs, f.coeffs[1])
            push!(rexps, f.exps[1])
            f = FreeAssAlgElem{T}(R, f.coeffs[2:end], f.exps[2:end], length(f) - 1)
        end
    end
    return rep_dict, FreeAssAlgElem{T}(R, rcoeffs, rexps, length(rcoeffs))
end

