using DataFrames
using CSV

export DataTable, add_data!, data_table, add_names!, add_name!, get_row, save, load_dt
export _key_to_name, _insert_or_overwrite


mutable struct DataTable
    path_list::Vector{String}
    data::DataFrame
    DataTable() = new(collect(_data_path_iterator()), DataFrame(:Name => String[]))
    DataTable(path_list::Vector{String}, data::DataFrame) = new(path_list, data)
end

data_table() = DataTable()
data_table(data::DataFrame) = DataTable(collect(_data_path_iterator()), data)


function add_data!(dt::DataTable, path::String)
    df = dt.data
    name = String(_name_from_savepath(path))
    dct  = load_dict(path)
  
    M = matroid_from_matroid_hex(name)
    add_name!(dt, name)

    row = get_row(df, name)

    for key in keys(dct)
      !occursin("Aut",key) && continue
        
        v = split(key, "_")
        if length(v) == 3 #Not a timed thing
          y = dct[key];
          I = ideal(y);
          I.gb = Oscar.IdealGens(y);
          

          df = _insert_or_overwrite(df, row, name, key, QuantumAutomorphismGroups.is_commutative(I))
          df = _insert_or_overwrite(df, row, name, key * "_max_deg", maximum(AbstractAlgebra.total_degree.(y)))
          continue
        end
        if length(v) == 4

          df = _insert_or_overwrite(df, row, name, key, dct[key])
          continue
        end
    end

    dt.data = df
    return df
end

function _insert_or_overwrite(df::DataFrame, row::DataFrameRow, name::String, key::String, value::Any)
  if !(key in names(df))
    df = outerjoin(df, DataFrame("Name" => name, key => value), on = :Name)
    return df
  end 
  row[key] = value
  return df
end

function add_data!(dt::DataTable, paths::Vector{String})
    for pth in paths
        add_data!(dt, pth)
    end
end
add_data!(dt::DataTable) = add_data!(dt, dt.path_list)

function get_row(df::DataFrame, name::String)
    return df[findfirst(==(name), df.Name), :]
end

get_row(dt::DataTable, name::String) = get_row(dt.data, name)

Base.getindex(dt::DataTable, name::String) = get_row(dt, name)
Base.getindex(dt::DataTable, M::Matroid) = get_row(dt, matroid_hex(M))

function Oscar.save(dt::DataTable, path::String="../data/data_table.csv")

    CSV.write(path, dt.data)
end

function load_dt(path::String="../data/data_table.csv")
    dt = DataTable()
    dt.data = CSV.read(path, DataFrame)
    return dt 
end

function add_names!(dt::DataTable)
    for pth in dt.path_list
      name = _name_from_savepath(pth)
      add_name!(dt, String(name))
    end

    return dt.data
end

function add_name!(dt::DataTable, name::String)
  df = dt.data
  row = filter(:Name => ==(name), df) 
  if size(row)[1] == 0
    M = matroid_from_matroid_hex(name)
    append!(df, DataFrame("Name" => name),promote = true,cols = :union)
  end
  
  return df
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
dt = data_table()
add_names!(dt)
add_data!(dt)

y = load_dict(uniform_matroid(1,3))["Aut_bases_-1"]
max(deg.(y)) rw


dt = load_dt();
for row in eachrow(dt.data)
  name = row.Name
  

end


df2 = select(dt.data, ["Name", "Aut_bases_-1", "Aut_circuits_-1","Aut_circuits_6", "Aut_bases_6"]);
df2 = filter(row -> row["Aut_bases_6"] !== missing || !ismissing(row["Aut_bases_-1"]), df2)
show(df2)
show(dt.data, allcols = true)
dt2 = data_table(df2)


r2 = dt2[fano_matroid()]
r3 = dt2[non_fano_matroid()]
select(DataFrame(r3), ["Name",r"bases"])

pth =save_path(matroid_hex("r0n1_1"))
pth2 = save_path(matroid_from_matroid_hex("r0n2_1"))

add_data!(dt, pth)
add_data!(dt, pth2)

name = _name_from_savepath(pth)
df[findfirst(==(name), df.Name), :]

"Name" in names(df)

name = String(_name_from_savepath(pth))
row = filter(:Name => ==(name), dt.data) 

=#

