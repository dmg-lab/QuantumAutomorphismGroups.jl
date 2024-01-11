
export computeGbOfMatroid,
    computeLpGbOfMatroid,
    getTiming

@doc raw"""

    computeGbOfMatroid(M::Matroid,structure::Symbol=:bases)

Compute the Groebner basis of the quantum automorphism group of a matroid for a given structure. The result is saved in the data folder. I will load the Gb from the data folder if it has already been computed.

# Examples
```julia
computeGbOfMatroid(uniform_matroid(3,4))

```
"""
function computeGbOfMatroid(M::Matroid,structure::Symbol=:bases)
    data_dir = "../data/"
    infoFile = ".info"

    name = getName(M) 
    path = "r$(rank(M))n$(length(M))/"* name
    folder, filename = split(path, "/")
    fullpath = data_dir * path * infoFile

    #Check if folder exists, if not create it
    isdir(data_dir * folder) || mkdir(data_dir * folder)


    #Check if info exists, if not create it
    if isfile(fullpath) 
        info = loadDict(fullpath)
    else
        info = Dict{String,Any}(
        "revlex_basis_encoding" =>String(M.pm_matroid.REVLEX_BASIS_ENCODING),
        )
    end
    if !haskey(info, "Aut_" * String(structure)) 
        #Compute
        gns , _ , _ = getMatroidRelations(M,structure)
        println("Computing Aut_$(String(structure)) for  $(name)") 
        gb, elapsed = @timed AbstractAlgebra.groebner_basis(gns)
        #Saving 
        info["Aut_" * String(structure)] = (gb);
        info["Aut_" * String(structure)*"_timed"] = elapsed
        saveDict(fullpath, info)

    else
        println("Already computed")
    end

    return info
end

function getTiming(M::Matroid,structure::Symbol=:bases)
    data_dir = "../data/"
    infoFile = ".info"

    name = getName(M) 
    path = "r$(rank(M))n$(length(M))/"* name
    folder, filename = split(path, "/")
    fullpath = data_dir * path * infoFile

    #Check if folder exists, if not create it
    isdir(data_dir * folder) || mkdir(data_dir * folder)

    #Check if info exists, if not compute it
    if isfile(fullpath) 
        info = loadDict(fullpath)
    else
        info = Dict{String,Any}(
        "revlex_basis_encoding" =>String(M.pm_matroid.REVLEX_BASIS_ENCODING),
        )
    end
    if !haskey(info, "Aut_" * String(structure)*"_timed") 
        if !haskey(info, "Aut_" * String(structure)) 
            info = computeGbOfMatroid(M,structure) 
            return info
        else
            println("Adding timing")
            #Compute
            gns , _ , _ = getMatroidRelations(M,structure)
            println("Computing Aut_$(String(structure)) for  $(name)") 
            _, elapsed = @timed AbstractAlgebra.groebner_basis(gns)
            #Saving 
            info["Aut_" * String(structure)*"_timed"] = elapsed
            saveDict(fullpath, info)
        end 
    else
        println("Already computed")
    end
end

getTiming(Name::String,structure::Symbol=:bases) = getTiming(nameToMatroid(Name),structure)
#=
computeGbOfMatroid(uniform_matroid(1,2),[:bases,:circuits])
=#
function computeGbOfMatroid(M::Matroid, strcts::Vector{Symbol})
    ans = Dict{String,Any}()
    for strct in strcts
        ans = computeGbOfMatroid(M,strct)
    end
    return ans
end


function toFreeAssAlgElem(U::Matrix{FreeAssAlgElem{QQFieldElem}}, p::Singular.slpalg{Singular.n_Q})
    coeffs = Oscar.QQ.(collect(Oscar.coefficients(p)))
    new_p = 0

    for (i, exp) in enumerate(collect(Oscar.Singular.exponent_words(p)))
        mon = coeffs[i]
        for i in exp
            mon *= U[i]
        end
        new_p += mon
    end
    return new_p
end

#=
computeLpGbOfMatroid(uniform_matroid(3,4),:bases)
=#
function computeLpGbOfMatroid(M::Matroid, structure::Symbol, n::Int=3)

    data_dir = "../data/"
    infoFile = ".info"

    name = getName(M) 
    path = "r$(rank(M))n$(length(M))/"* name
    folder, filename = split(path, "/")
    fullpath = data_dir * path * infoFile

    #Check if folder exists, if not create it
    isdir(data_dir * folder) || mkdir(data_dir * folder)


    #Check if info exists, if not create it
    if isfile(fullpath) 
        info = loadDict(fullpath)
    else
        info = Dict{String,Any}(
        "revlex_basis_encoding" =>String(M.pm_matroid.REVLEX_BASIS_ENCODING),
        )
    end
    if !haskey(info, "Aut_" * String(structure) * "_fast") 
        #Compute
        println("Computing Aut_$(String(structure)) for  $(name)") 
        gns, U , A = getMatroidRelations(M,structure)
        I = Oscar.ideal(A,gns)
        Oscar.groebner_assure(I,n);
        Oscar.singular_assure(I.gb,n)
        s_gb = gens(I.gb.S)
        gb =  AbstractAlgebra.groebner_basis(map(x->toFreeAssAlgElem(U,x),s_gb))

        #Saving 
        info["Aut_" * String(structure)*"_fast"] = (gb);
        saveDict(fullpath, info)

    else
        println("Already computed")
    end

    return info

end













