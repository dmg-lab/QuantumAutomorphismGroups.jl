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
        #Saving 
        info["Aut_" * String(structure)] = (gb);
        saveDict(fullpath, info)

    else
        println("Already computed")
    end

    return info
end

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
    if !haskey(info, "Aut_" * String(structure)) 
        #Compute
        println("Computing Aut_$(String(structure)) for  $(name)") 
        gns, U , A = getMatroidRelations(M,structure)
        I = Oscar.ideal(A,gns)
        Oscar.groebner_assure(I,n);
        Oscar.singular_assure(I.gb,n)
        s_gb = gens(I.gb.S)
        gb =  AbstractAlgebra.groebner_basis(map(x->toFreeAssAlgElem(U,x),s_gb))

        #Saving 
        info["Aut_" * String(structure)] = (gb);
        saveDict(fullpath, info)

    else
        println("Already computed")
    end

    return info

end



#=






#Temporary computation, include what should be computated
M = uniform_matroid(3,5)
computeLpGbOfMatroid(M,:bases)
=#

##= Database computation


global Droids = []
for n in 1:7, r in 1:n

    db = Polymake.Polydb.get_db()
    collection = db["Matroids.Small"]

    cursor=Polymake.Polydb.find(collection, Dict("RANK" => r,"N_ELEMENTS"=>n))

    append!(Droids,Matroid.(cursor))



     
end
sort!(Droids,by=x->length(getMatroidRelations(x,:bases)[1]))


for M in Droids
    computeLpGbOfMatroid(M,:bases)
end


#=#


#=



db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]

cursor=Polymake.Polydb.find(collection, Dict("RANK"=>3, "N_ELEMENTS"=>6))
Droids = Matroid.(cursor)

sort!(Droids,by=x->length(getMatroidRelations(x,:bases)[1]))





strcts = Symbol[:bases,:circuits]

for M in Droids
    computeLpGbOfMatroid(M,:bases)
    println("Done with $(getName(M))")
end

info = computeGbOfMatroid(Droids[10])
I = info["Aut_bases"]
isCommutative(I)
I

computeGbOfMatroid(uniform_matroid(3,4),:circuits)
Oscar.ideal(A,gns)

M = uniform_matroid(3,4)
gns, U , A = getMatroidRelations(M,:bases)


=#





#=


collect(Oscar.Singular.exponent_words(gb[23]))

methodswith(typeof(gb[1]))


FreeAssAlgElem(gb[1])
# convert to FreeAssAlgebraElem?


I = Oscar.ideal(A,gns);
ideal_membership(gns[1],I,3)


gens(I)



=#







