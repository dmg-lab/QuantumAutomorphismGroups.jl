using Oscar
import AbstractAlgebra: Generic.AhoCorasickAutomaton


function toNamedTuple(v::Vector{Dict{Int,Int}})
    res =Vector{NamedTuple}()
    for ds in v 
        tmp = (; (Symbol(k) => v for (k,v) in ds)...)
        push!(res, NamedTuple(tmp))
    end
    return Tuple(res)
end

function toDict(v::Vector{NamedTuple})
    res =Vector{Dict{Int,Int}}()
    for ds in v 
        dct = Dict{Int,Int}()
        for (k,v) in pairs(ds)
            dct[parse(Int,string(k))] = Int(v)
        end
        push!(res, dct)
    end
    return res
end




function saveAhoCorasick(path::String,Aho::AhoCorasickAutomaton)
   gt = toNamedTuple(Aho.goto)   
   save(path,(gt,Aho.fail,Aho.output))
end

function loadAhoCorasick(path::String)
    res = load(path)
    l = toDict(collect(res[1]))
    return AhoCorasickAutomaton(l, res[2], res[3])
end

function saveDict(path::String, dct::Dict)
    tmp = (; (Symbol(k) => v for (k,v) in dct)...)
    save(path, NamedTuple(tmp))
end

function loadDict(path::String)
    namedTuple = load(path)
    dct = Dict{String,Any}()
    for (k,v) in pairs(namedTuple)
        dct[string(k)] = v
    end
    return dct
end



#=
M = uniform_matroid(2,3)

info = Dict(
"revlex_basis_encoding" =>String(M.pm_matroid.REVLEX_BASIS_ENCODING),
"hasSymmetries" => true
)

saveDict("./uM.info",info)
loadDict("./uM.info") 
=#
