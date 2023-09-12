include("./quantumMatroid.jl")


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
        gb = AbstractAlgebra.groebner_basis(gns)
        println("Computed")
        #Saving 
        info["Aut_" * String(structure)] = (gb);
        saveDict(fullpath, info)

    else
        println("Already computed")
    end

    return info
end

function computeGbOfMatroid(M::Matroid, strcts::Vector{Symbol}=Symbol[:bases,:flats,:circuits,:rank])
    ans = Dict{String,Any}()
    for strct in strcts
        ans = computeGbOfMatroid(M,strct)
    end
    return ans
end



#=

#Temporary computation, include what should be computated
M = uniform_matroid(0,2)
computeGbOfMatroid(M,Symbol[:bases,:flats,:circuits,:rank])


info = computeGbOfMatroid(uniform_matroid(2,5),:circuits)


=#

#= Database computation
Droids = []
for n in 1:5, r in 1:n

    db = Polymake.Polydb.get_db()
    collection = db["Matroids.Small"]

    cursor=Polymake.Polydb.find(collection, Dict("RANK" => r,"N_ELEMENTS"=>n))
    Droids = vcat(Droids,Matroid.(cursor))

end

Droids
strcts = Symbol[:bases,:flats,:circuits,:rank]

for M in Droids
    computeGbOfMatroid(M,:bases)
end


=#







