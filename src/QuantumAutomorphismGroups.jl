
module QuantumAutomorphismGroups
    using Oscar

    DATA_DIR = "../data/"
    INFO_FILETYPE = ".info"


    include("./is_commutative.jl")
    include("./compute_and_store_gb.jl")
    include("./data_table.jl")
    include("./handle_old_data.jl")
end


