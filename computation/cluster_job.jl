using Oscar
using QuantumAutomorphismGroups

using CSV
using DataFrames

csv="../data/data_table.csv"

df = CSV.read(csv, DataFrame) 
for name in df.Name
    M = nameToMatroid(String(name))
    if length(M) >= 5
        continue
    end
    computeLpGbOfMatroid(M,:bases)
    println("Done with $name")
end 

#=
global Droids = []
for n in 7:7, r in 1:n

    db = Polymake.Polydb.get_db()
    collection = db["Matroids.Small"]

    cursor=Polymake.Polydb.find(collection, Dict("RANK" => r,"N_ELEMENTS"=>n))

    append!(Droids,Matroid.(cursor))
end

sort!(Droids,by=x->length(getMatroidRelations(x,:bases)[1]))


for M in Droids
    computeGbOfMatroid(M,:bases)
end
=#
#=
R = "0000*******************************"
M = matroid_from_revlex_basis_encoding(R,3,7)
N = matroid_from_nonbases(nonbases(M),6)
loadInfo(N)

M=uniform_matroid(1,3)
computeGbOfMatroid(M,:bases)

computeGbOfMatroid(M,:bases)

=#



