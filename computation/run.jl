Pkg.activate("../")

using Oscar
include("./utils/save.jl")




#Temporary
matroidName = "non_fano_matroid"
M = non_fano_matroid()
structure = :bases



data_dir = "./data/"
infoFile = "stats.info"

M = @isdefined(M) ? M : uniform_matroid(1,2) 


path = "r$(rank(M))n$(length(M))/"* matroidName



folder, filename = split(path, "/")
fullpath = data_dir * folder * "/" * filename * "/"

#Check if folder exists, if not create it
isdir(data_dir * folder) || mkdir(data_dir * folder)
isdir(fullpath) || mkdir(fullpath)


#Check if info exists, if not create it
if isfile(fullpath * infoFile ) 
    info = loadDict(fullpath * infoFile)
else
    info = Dict{String,Any}(
    "revlex_basis_encoding" =>String(M.pm_matroid.REVLEX_BASIS_ENCODING),
    )
end


#Temporary computation, include what should be computated
if !haskey(info, "Aut_" * String(structure)) 
    #Compute
    include("./utils/getData.jl");
    #Saving 
    info["Aut_" * String(structure)] = (gens, dur);
    println("Saving Infos")
    saveDict(fullpath * infoFile, info)
    #Save the ahocorasick automaton
    println("Saving AhoCorasick")
    saveAhoCorasick(fullpath*"Aut_"* String(structure) * ".aho", aut)
else
    println("Already computed")
end
    







