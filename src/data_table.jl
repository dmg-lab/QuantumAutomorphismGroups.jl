using DataFrames
using CSV

export DataTable, add_data!
export _key_to_name


struct DataTable
    path_list::Vector{String}
    data::DataFrame
    DataTable() = new(collect(_data_path_iterator()), DataFrame())
end

data_table = DataTable()


function add_data!(df::DataFrame, path::String)
    name = _name_from_savepath(path)
    dct  = load_dict(path)

    old_row = filter(:Name => ==(name), df)
    
    if size(old_row,1) == 0
        data = Dict{String,Any}("Name" => name)
    else
        data = Dict{String,Any}(map(x -> String(x[1]) => x[2], pairs(old_row[1,:])))
    end

    for keys in keys(dct)
        data[keys] = dct[keys]
    end
end


function _key_to_name(key::String)
    word_vector = split(key, "_")
    @assert word_vector[1] == "Aut" "The key does not start with Aut"
    if length(word_vector) == 3
      return uppercasefirst(word_vector[2]) * "_" * word_vector[3]
    end
    return uppercasefirst(word_vector[2]) * "_" * word_vector[3] * "_" * word_vector[4]
end

#=
using Oscar
_data_path_iterator()

for path in _data_path_iterator()
  dct = load_dict(path)
  if haskey(dct, "Aut_bases_6")
    y = dct["Aut_bases_6"];
    I = ideal(y);
    I.gb = Oscar.IdealGens(y);
    println(QuantumAutomorphismGroups.is_commutative(I))
  end

end

M = fano_matroid()
dct = load_dict(M)

is_commutative(dct["Aut_bases_6"])

QuantumAutomorphismGroups.is_commutative(dct["Aut_bases_6"])
#captilize the first letter of the key "matroid"
str = "matroid"
capitalized_str = uppercasefirst(str)


using DataFrames


# Example usage
dt = DataFrame()
_name_from_save_path(save_path(M))


dct = Dict("matroid" => "hex", "rank" => 3, "n" => 6)

new_row = DataFrame(dct)
add!

add_data!(dt, dct)



=#
