using Oscar
using QuantumAutomorphismGroups

using CSV
using DataFrames
#=
#M = non_fano_matroid()
#computeLpGbOfMatroid(M,:bases)
task_number = parse(Int,ENV["SLURM_ARRAY_TASK_ID"])

n_tasks = parse(Int,ENV["SLURM_ARRAY_TASK_COUNT"])
println("Task number: ", task_number)


global Droids = []
for n in 7:7
    db = Polymake.Polydb.get_db()
    collection = db["Matroids.Small"]
    cursor=Polymake.Polydb.find(collection, Dict("N_ELEMENTS"=>n))
    append!(Droids,Matroid.(cursor))
end


global myDroids = view(Droids,task_number:n_tasks:length(Droids))
    

for M in myDroids
    computeLpGbOfMatroid(M,:bases)
    println("Computed for $(getName(M))")
end

=#

M = fano_matroid()
computeLpGbOfMatroid(M,:bases)



