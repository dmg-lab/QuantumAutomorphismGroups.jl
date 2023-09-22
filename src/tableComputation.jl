using DataFrames
using CSV

include("quantumMatroid.jl")

#=

df = addAutToDF()
names(df)
=#
function addAutToDF(recompute::Bool=false)

    dt = loadAll()
    if isfile("../data/data_table.csv")
        df = CSV.read("../data/data_table.csv",DataFrame)
    else
        df = DataFrame(Name=String[])
    end
    
    for (name, dict) in dt
        
        row = filter(:Name => ==(name),df)
        size(row)[1] == 0 ? 
        (data = Dict{String,Any}(); data["Name"] = name) : 
        data = Dict{String,Any}(map(x->String(x[1])=>x[2],pairs(row[1,:])))

        if !haskey(data,"rank") ||  !haskey(data,"length") || recompute 
            println("Computing rank and length of $name")
            M = nameToMatroid(name)
            data["rank"] = rank(M)
            data["length"] = length(M)
        end


        
        for key in keys(dict)
            #if starts with Aut, then add to df
            if occursin("Aut",key)
                (_,structure) = split(key,"_")
                if !recompute && haskey(data,"Aut_$(uppercase(structure[1]))")
                    continue
                end
                data["Aut_$(uppercase(structure[1]))"] = isCommutative(dict[key])[1]
                data["Aut_$(uppercase(structure[1]))_stored"] = true
            end 
        end
        #only appand if this name does not exist yet
        if size(row)[1] == 0
            append!(df,data,promote=true,cols=:union)
            continue
        end
        new_line = select!(DataFrame(data),Symbol.(names(df)))

        df[df.:Name .== name, :] = new_line


    end


    CSV.write("../data/data_table.csv",df)


    return df    
end

#=
df = DataFrame(Name=String[],Simple=Bool[])
addToDf(df,uniform_matroid(1,3),"test","random")
=#
function addToDf(df::DataFrame, M::Matroid, dataPoint::Any, colname::String, recompute::Bool=false)
    name = getName(M)

    row = filter(:Name => ==(name),df)
    size(row)[1] == 0 ? 
    (data = Dict{String,Any}(); data["Name"] = name) : 
    data = Dict{String,Any}(map(x->String(x[1])=>x[2],pairs(row[1,:])))


    if !recompute && haskey(data,colname)
        return df
    end

    data[colname] = dataPoint
    append!(df,data,promote=true,cols=:union)
    return df
end





#=
loadInfo(uniform_matroid(1,3))

loadInfo("r1n2_1")

=#

function loadInfo(name::String)

    data_dir = "../data/"
    infoFile = ".info"

    folder = split(name,"_")[1]
    
    path = folder * "/"* name
    fullpath = data_dir * path * infoFile

    #Check if folder exists, if not create it
    isdir(data_dir * folder) || throw(ArgumentError("No data for this matroid exists"))


    #Check if info exists, if not create it
    if isfile(fullpath) 
         return  loadDict(fullpath)
    else
        throw(ArgumentError("No data for this matroid exists"))
    end

end


loadInfo(M::Matroid) = loadInfo(getName(M))


#=
loadAll()
=#
function loadAll()
    data_dir = "../data/"
    infoFile = ".info"

    data = Dict{String,Dict{String,Any}}()
    
    for folder in readdir(data_dir)
        if isdir(data_dir * folder)
            for file in readdir(data_dir * folder)
                if isfile(data_dir * folder * "/" * file)
                    #cut of the .info
                    fileName = split(file,".")[1]
                    data[fileName] = loadDict(data_dir * folder * "/" * file)
                end
            end
        end
    end
    return data
end



function toDataFrame(D::Dict{String,Dict{String,Bool}})
    otp = DataFrame()

    for key in keys(Test)
        df = DataFrame(Test[key])
        df[!,:Name] .= key
        df = select!(df,[:Name,Symbol.(keys(Test[key]))...])
        append!(otp,df,promote=true,cols=:union)
    end
    return otp
end








#= Save and Load File

path = "../data/data_table.csv"
df = CSV.read(path,DataFrame)
select(df,[:Name,:length,:rank, :Aut_B,:Aut_R,:Aut_F,:Aut_C,:describe])

=#
