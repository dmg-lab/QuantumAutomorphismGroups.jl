#To get the old data

export old_save_path, old_data_path_iterator, compare_new_and_old



function old_save_path(M::Matroid)
  return "$(DATA_DIR)r$(rank(M))n$(length(M))/$(matroid_hex(M))$(INFO_FILETYPE)"
end


function old_data_path_iterator()
    DATA_DIR = "../data_old/"
    folders = readdir(DATA_DIR)
    return (DATA_DIR * folder * "/" * file for folder in folders if isdir(DATA_DIR * folder) for file in readdir(DATA_DIR * folder))
end

function compare_new_and_old()
  for old_path in old_data_path_iterator()
    old_dict = load_dict(old_path)
    name = _name_from_savepath(old_path)
    save_path = QuantumAutomorphismGroups.save_path(matroid_from_matroid_hex(name))
    
    #check if save_path exists
    if !isfile(save_path) 
      println("Creating new file $save_path")
      new_dict = Dict{String,Any}()
    else
      new_dict = load_dict(save_path)
    end

    
    for (key,val) in pairs(old_dict)

      !occursin("Aut",key) && continue
        
      v = split(key, "_")

      if length(v) == 2
        structure = v[2]
        name_in_new_dct = "Aut_" * structure * "_-1"
        timing_name = "Aut_" * structure *  "_-1" *  "_timed"
        if !haskey(new_dict, name_in_new_dct)
          println("Adding $name_in_new_dct to $save_path")
            new_dict[name_in_new_dct] = val
            if haskey(old_dict, "Aut_" * structure * "_timed")
              new_dict[timing_name] = old_dict["Aut_" * structure * "_timed"]
            end
        end
      end
    end
  save_dict(save_path, new_dict)



  end


end





#=
compare_new_and_old()


dt = data_table()
m = matroid_from_matroid_hex("r3n6_00001")
load_dict(m)


dt.path_list = collect(old_data_path_iterator())

load_dict(dt.path_list[45])

add_names!(dt)
add_data!(dt)

dt_new = load_dt()
=#


