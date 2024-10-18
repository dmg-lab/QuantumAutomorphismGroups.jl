export compute_and_store_gb
export save_dict, load_dict, save_path

@doc raw"""

    compute_and_store_gb(M::Matroid; structure::Symbol=:bases, deg_bound::Int = -1)

Compute the Groebner basis of the quantum automorphism group of a matroid for a given structure. The result is saved in the data folder. I will load the Gb from the data folder if it has already been computed.

# Examples
```julia
computeGbOfMatroid(uniform_matroid(3,4))

```
"""
function compute_and_store_gb(M::Matroid, structure::Symbol=:bases; deg_bound::Int = -1)
    DATA_DIR = "../data/"
    INFO_FILETYPE = ".info"

    name = matroid_hex(M) 
    path = "r$(rank(M))n$(length(M))/"* name
    folder, filename = split(path, "/")
    fullpath = DATA_DIR * path * INFO_FILETYPE

    #Check if folder exists, if not create it
    isdir(DATA_DIR * folder) || mkdir(DATA_DIR * folder)
    
    computation_name = "Aut_" * String(structure) * "_" * string(deg_bound)



    #Check if info exists, if not create it
    if isfile(fullpath) 
        info = load_dict(fullpath)
    else
        info = Dict{String,Any}(
        )
    end
    if !haskey(info, computation_name) || !haskey(info, computation_name * "_timed")
        #Compute
        I = quantum_automorphism_group(M, structure)
        println("Computing Aut_$(String(structure)) for $(name) with degree bound $(deg_bound)")

        gb, elapsed = @timed groebner_basis(I, deg_bound)
        #Saving 
        info[computation_name] = collect(gb);
        info[computation_name*"_timed"] = elapsed
        save_dict(fullpath, info)

    else
        println("Already computed")
    end

    return info
end


#=

All of these are functions that are used to save and load the groebner basis information in a dictionary format.


using Oscar
M = uniform_matroid(1,3)
save_path(M)

load_dict(M)
=#
function to_named_tuple(v::Vector{Dict{Int,Int}})
    res =Vector{NamedTuple}()
    for ds in v 
        tmp = (; (Symbol(k) => v for (k,v) in ds)...)
        push!(res, NamedTuple(tmp))
    end
    return Tuple(res)
end

function to_dict(v::Vector{NamedTuple})
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

function save_dict(path::String, dct::Dict)
    tmp = (; (Symbol(k) => v for (k,v) in dct)...)
    save(path, NamedTuple(tmp))
end

function load_dict(path::String)
    namedTuple = load(path)
    dct = Dict{String,Any}()
    for (k,v) in pairs(namedTuple)
        dct[string(k)] = v
    end
    return dct
end

function load_dict(M::Matroid)
  return load_dict(save_path(M))
end

function save_path(M::Matroid)
  return "$(DATA_DIR)r$(rank(M))n$(length(M))/$(matroid_hex(M))$(INFO_FILETYPE)"
end

