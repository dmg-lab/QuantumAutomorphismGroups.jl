using DataFrames
using CSV

include("quantumMatroid.jl")


function addAutToDF(df::DataFrame, M::Matroid, structure::Symbol=:bases, alt::Bool=false; n::Int=3, recompute::Bool=false)
    name = getName(M)

    #check if a column with the name "Aut_$(String(structure))" already exists
    colname = "Aut_$(uppercase(String(structure)[1]))"
    data = Dict{String,Any}()
    data["Name"] = name

    row = filter(:Name => ==(name),df)
    size(row)[1] == 0 ? 
    (data = Dict{String,Any}(); data["Name"] = name) : 
    data = Dict{String,Any}(map(x->String(x[1])=>x[2],pairs(row[1,:])))

    if !recompute && haskey(data,colname)
        return df
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

addAutToDF(df,uniform_matroid(1,3),:circuits)





true ?
(data = Dict{String,Any}();data["Name"]="test") :
data = Dict{String,Any}(df[1,:])
data
=#

#= Save and Load File

path = "../data/table.csv"
CSV.write(path,df)


new_df = CSV.read(path,DataFrame)
=#

