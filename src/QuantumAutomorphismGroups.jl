
module QuantumAutomorphismGroups
    using Oscar

    DATA_DIR = "../data/"
    INFO_FILETYPE = ".info"


    include("./is_commutative.jl")
    include("./compute_and_store_gb.jl")
end


