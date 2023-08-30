Pkg.activate("../")


using Oscar



function splitString(s::String, n::Int)
    if length(s) <= n
        return [s]
    else
        return [s[1:n] ; splitString(s[n+1:end], n)]
    end
end

toHex(s::String) = string(parse(Int, s, base=2), base=16)

function getMatroidId(M::Matroid)
    rev = replace(String(M.pm_matroid.REVLEX_BASIS_ENCODING), "*"=>"1")
    sRev = splitString(rev, 8)
    sHex = map(toHex, sRev)
    return reduce(*, sHex)
end  

#=
M = uniform_matroid(5, 10) 
getMatroidId(M)
=#




