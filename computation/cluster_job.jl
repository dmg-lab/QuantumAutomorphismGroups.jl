using Oscar
using QuantumAutomorphismGroups

using CSV
using DataFrames

#=
M = uniform_matroid(1,3)
compute_and_store_gb(M,:bases,deg_bound=6)
load_dict(M)
=#


task_number = parse(Int,ENV["SLURM_ARRAY_TASK_ID"])

n_tasks = parse(Int,ENV["SLURM_ARRAY_TASK_COUNT"])
println("Task number: ", task_number)


global Droids = []
for n in 1:5
    db = Polymake.Polydb.get_db()
    collection = db["Matroids.Small"]
    cursor=Polymake.Polydb.find(collection, Dict("N_ELEMENTS"=>n))
    append!(Droids, Matroid.(cursor))
end


global myDroids = view(Droids,task_number:n_tasks:length(Droids))
    
#bases, circuits, flats
for M in myDroids
    dg_bound = -1
    compute_and_store_gb(M, :bases; deg_bound=dg_bound)
    println("Computed for $(matroid_hex(M)) with degree bound $(dg_bound)")
end





