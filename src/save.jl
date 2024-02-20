using Oscar
using Oscar.JSON
using Oscar.AbstractAlgebra
using Oscar.AbstractAlgebra.Generic

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

#=
function strToAho(str::String)
    x = JSON.parse(str)
    newGoto = Vector{Dict{Int,Int}}()
    for d in  x["goto"]
        newdct = Dict{Int,Int}() 
        for (k,v) in d
            newdct[parse(Int,k)] = v
        end
        push!(newGoto,newdct)
    end
    newFail = Vector{Int}(x["fail"])
    newOtp = Vector{Tuple{Int,Vector{Int}}}()
    for ele in x["output"]
        push!(newOtp,(ele[1],Vector{Int}(ele[2])))
    end
    return  AbstractAlgebra.Generic.AhoCorasickAutomaton(newGoto,newFail,newOtp)
end

function saveAhoCorasick(path::String,aho::AbstractAlgebra.Generic.AhoCorasickAutomaton)
    str = JSON.json(aho)
    open(path, "w") do file
        write(file, str);
    end;
    return
end
function loadAhoCorasick(path::String)
    str = read(path,String)
    return strToAho(str)
end

=#


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
