using DataFrames
using CSV

include("quantumMatroid.jl")


function addAutToDF(df::DataFrame, M::Matroid, structure::Symbol=:bases, alt::Bool=false; n::Int=3, recompute::Bool=false)
    name = getName(M)

    #check if a column with the name "Aut_$(String(structure))" already exists
    colname = "Aut_$(uppercase(String(structure)[1]))"
    data = Dict{String,Any}()
    data["Name"] = name
    #check if a row with this name already exists
    for i in 1:size(df)[1]
        if df[i,:Name] == name
            #add column if it does not exist 
            filter(x->x==colname, names(df)) == [] ? df[!,Symbol(colname)] = false : nothing




            if recompute
                df[i,Symbol(colname)] = isCommutative(M,structure,alt,n)[1]
                return df
            else
                return df
            end
        end
    end
    data[colname] = isCommutative(M,structure,alt,n)[1]
    append!(df,data,promote=true,cols=:union)
    return df
end


function addToDf(df::DataFrame, M::Matroid, data::Any, colname::Symbol)
    #check if a row with this name already exists
    for i in 1:size(df)[1]
        if df[i,:Name] == name
             
            df[i,colname] = data
            return df
        end
    end





end






function getNewRow(df::DataFrame, M::Matroid, alt::Bool=false; n::Int=3)
    aut_b = isCommutative(M,:bases,alt,n)[1]
    name = getName(M)
    #check if a row with this name already exists
    for i in 1:size(df)[1]
        if df[i,:Name] == name
            df[i,:AutB] = aut_b
            return df
        end
    end
    push!(df,(getName(M),isCommutative(M,:bases,alt,n)[1]))
    return df
end


#=
df = DataFrame(Name=String[],Simple=Bool[])

addAutToDF(df,uniform_matroid(3,4),:bases)
=#

#= Save and Load File

path = "../data/table.csv"
CSV.write(path,df)


new_df = CSV.read(path,DataFrame)
=#

