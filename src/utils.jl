using Oscar

export getName,
    nameToRevlex,
    nameToMatroid,
    powerset,
    custom_deg,
    altName




function matroidIncidence(M::Matroid,t::Symbol=:circuits)
    rows = eval(t)(M);

    groundset = matroid_groundset(M);
    incidence = zeros(Int64,length(rows),length(groundset));

    for i in 1:length(rows)
        for j in 1:length(groundset)
            if groundset[j] in rows[i]
                incidence[i,j] = 1;
            end
        end
    end
    
    return incidence;
end


function toAdjacency(Inc::Matrix{Int64})
    t = hcat(zeros(Int64,size(Inc,1),size(Inc,1)),Inc)
    t1 = hcat( transpose(Inc),zeros(Int64,size(Inc,2),size(Inc,2)))
    adj = vcat(t,t1)
    return adj
end


function Base.contains(V::Vector,v::Vector)
    isempty(v)  && return false
    S = string(V)
    s = string(v) 
    x = findfirst("[",s)
    s = s[x[1]+1:end-1]
    return contains(S,s)
end


function getFreePercentageOfMemory()
    free_mem = Sys.free_memory() / 2^20
    total_mem = Sys.total_memory() / 2^20
    return free_mem/total_mem
end

function memoryCritical()
    return getFreePercentageOfMemory() < 0.15
end



function deleteContaining!(V::Vector,verbose::Bool=false)
     
    isempty(V) && return V
    for size in 1:length(V)
        i = 1 
        while i <= length(V[size])
            for bigger in size+1:length(V)
                isempty(V[bigger]) && continue
                j = 1
                while j <= length(V[bigger])
                    if  contains(V[bigger][j],V[size][i]) 
                        splice!(V[bigger],j) 
                    else
                        j+=1
                    end
                end
            end
        i+=1
        end
    end
    return V
end


function splitString(s::String, n::Int)
    if length(s) <= n
        return [s]
    else
        return [s[1:n] ; splitString(s[n+1:end], n)]
    end
end

toHex(s::String) = string(parse(Int, s, base=2), base=16)
toBin(s::String) = string(parse(Int, s, base=16), base=2)

function getName(M::Matroid)
    rev = replace(String(M.pm_matroid.REVLEX_BASIS_ENCODING), "*"=>"1")
    sRev = splitString(rev, 8)
    for i in 1:length(sRev)
        hex_str = toHex(sRev[i])
        if i < length(sRev) && length(hex_str) < 2
            hex_str = "0"^(2-length(hex_str))*hex_str
        end
        sRev[i] = hex_str
    end
    return "r$(rank(M))n$(length(matroid_groundset(M)))"*"_"*reduce(*, sRev)
end  

popfirst!(s::String) = s[2:end]

function altName(M::Matroid)
    rev = replace(String(M.pm_matroid.REVLEX_BASIS_ENCODING), "*"=>"1")
    strng = toHex(rev)
    supposed_length = Int(ceil(binomial(length(M),rank(M)) / 4 ))
    if supposed_length < length(strng)
        return strng
    else
        hex_str = "0"^(supposed_length-length(strng))*strng
        return hex_str
    end
end


#=
length(M)
altName(uniform_matroid(1,3))

Int(ceil(binomial(7,3)/4))
name = "r3n7_00c360dd7"
M = nameToMatroid(name)
altName(M)
M  = nameToMatroid("r2n6_000f")
rev = replace(String(M.pm_matroid.REVLEX_BASIS_ENCODING), "*"=>"1")
QuantumAutomorphismGroups.toHex(rev)
sRev = QuantumAutomorphismGroups.splitString(rev,8)
getName(M)
=#
function nameToRevlex(S::String, r::Int, n::Int)
    s = splitString(S, 2)
    for i in 1:length(s)
        bin_str = toBin(s[i])
        if i < length(s) && length(bin_str) < 8
            bin_str = "0"^(8-length(bin_str))*bin_str
        end 
        s[i] = bin_str
    end
    if (length(s)-1)*8 + length(s[end]) < binomial(n,r)
        s[end] = "0"^(binomial(n,r) - ((length(s)-1)*8 + length(s[end])))*s[end]
    end

    sBin = map(x->replace(x, "1"=>"*"), s)
    return reduce(*, sBin)
end    


function nameToMatroid(S::String)
    sep = split(S,"_")
    (r,n) = parse.(Int,split(sep[1][2:end],"n"))


    revl = nameToRevlex(String(sep[2]),r,n)
    matroid = matroid_from_revlex_basis_encoding(revl,r,n)
    return matroid
end    

function powerset(v::Vector{<:Integer}, min::Int =0, max::Int=length(v))
    ans = Vector{Int}[]
    for k in min:max
        push!(ans,Oscar.subsets(v,k)...)
    end
    return ans
end


custom_deg(p::FreeAssAlgElem{T}) where T<:FieldElem = maximum([length(x) for x in collect(exponent_words(p))])
custom_deg(p::Vector{<:FreeAssAlgElem{T}}) where T<:FieldElem = maximum([custom_deg(x) for x in p])

#=
M = fano_matroid() 
N = nameToMatroid(getName(M))
N.pm_matroid.REVLEX_BASIS_ENCODING == M.pm_matroid.REVLEX_BASIS_ENCODING

H = nameToMatroid("r1n2_1")
H.pm_matroid.REVLEX_BASIS_ENCODING

M.pm_matroid.REVLEX_BASIS_ENCODING == N.pm_matroid.REVLEX_BASIS_ENCODING
rank(M) == rank(N)
length(matroid_groundset(M)) == length(matroid_groundset(N))

getName(uniform_matroid(3,9))
=#




#function normal_form_with_rep(
#    f::FreeAssAlgElem{T},
#    g::Vector{FreeAssAlgElem{T}},
#    aut::AbstractAlgebra.Generic.AhoCorasickAutomaton,
#) where {T}
#rep_dict = Dict{Int, FreeAssAlgElem{T}}()
#    R = parent(f)
#    rexps = Monomial[]
#    rcoeffs = T[]
#    while length(f) > 0
#        ok, left, right, match_index = AbstractAlgebra.Generic.gb_divides_leftmost(f.exps[1], aut)
#        if ok
#            qi = AbstractAlgebra.divexact_right(f.coeffs[1], g[match_index].coeffs[1])
#            f = AbstractAlgebra.Generic._sub_rest(f, AbstractAlgebra.Generic.mul_term(qi, left, g[match_index], right), 1)
#            rep_dict[match_index] = AbstractAlgebra.Generic.mul_term(qi, left, g[match_index], right)
#        else
#            push!(rcoeffs, f.coeffs[1])
#            push!(rexps, f.exps[1])
#            f = FreeAssAlgElem{T}(R, f.coeffs[2:end], f.exps[2:end], length(f) - 1)
#        end
#    end
#    return rep_dict, FreeAssAlgElem{T}(R, rcoeffs, rexps, length(rcoeffs))
#end






