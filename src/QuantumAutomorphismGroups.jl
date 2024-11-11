
module QuantumAutomorphismGroups
    using Oscar

    DATA_DIR = "../data/"
    INFO_FILETYPE = ".info"


    include("./is_commutative.jl")
    include("./compute_and_store_gb.jl")
    include("./data_table.jl")

    #Not sure if this is the best place for this
    include("./handle_old_data.jl")
    include("./create_tex_table.jl")
end


