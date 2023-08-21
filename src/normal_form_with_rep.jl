
function normal_form_with_rep(
    f::FreeAssAlgElem{T},
    g::Vector{FreeAssAlgElem{T}},
    aut::AhoCorasickAutomaton,
) where {T}
rep_dict = Dict{Int, FreeAssAlgElem{T}}()
    R = parent(f)
    rexps = Monomial[]
    rcoeffs = T[]
    while length(f) > 0
        ok, left, right, match_index = gb_divides_leftmost(f.exps[1], aut)
        if ok
            qi = divexact(f.coeffs[1], g[match_index].coeffs[1])
            f = _sub_rest(f, mul_term(qi, left, g[match_index], right), 1)
            rep_dict[match_index] = mul_term(qi, left, g[match_index], right)
        else
            push!(rcoeffs, f.coeffs[1])
            push!(rexps, f.exps[1])
            f = FreeAssAlgElem{T}(R, f.coeffs[2:end], f.exps[2:end], length(f) - 1)
        end
    end
    return rep_dict, FreeAssAlgElem{T}(R, rcoeffs, rexps, length(rcoeffs))
end

