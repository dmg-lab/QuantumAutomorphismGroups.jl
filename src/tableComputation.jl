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

dt = loadAll()
row = String[]
inclusionMatrix = Vector[]
for (name,dict) in dt
    for key in keys(dict)
        if occursin("Aut",key)
            (_,structure) = split(key,"_")
            structure = uppercase(structure[1])
            push!(row,name * "_" * structure)
            v = Bool[] 
            for (name2,dict2) in dt
                for key2 in keys(dict2)
                    if occursin("Aut",key2)
                        #Check if the ideal behind key1 is contained in the ideal behind key2
                        isInIdeal(dict[key],dict2[key2]) ? push!(v,true) : push!(v,false)
                    end
                end

            end
            push!(inclusionMatrix,v)
        end
    end
end

save("./tst.tst",inclusionMatrix)


#= Save and Load File

path = "../data/data_table.csv"
df = CSV.read(path,DataFrame)
select(df,[:Name,:length,:rank, :Aut_B,:Aut_R,:Aut_F,:Aut_C,:describe])


df[!,"describe"] = Vector{Union{Missing,String}}(missing,size(df)[1])

for row in eachrow(df)
    name = String(row["Name"])
    M = nameToMatroid(name)
    rank(M)<1 && continue 
    row["describe"] = Oscar.describe(automorphism_group(M))
end
select(df,[:Name,:length,:rank,:describe])

df 


CSV.write(path,df)

M = uniform_matroid(0,10)
rank(M)

rank(M)<1 ? automorphism_group(dual_matroid(M)) : automorphism_group(M)

matroid_from_bases([[]],0)





x = automorphism_group(uniform_matroid(2, 2))
typeof(x)
symmetric_group(2)== automorphism_group(uniform_matroid(2, 2))



sort(df,[:length,:rank])
throw(ArgumentError("The Matroid has to be on at least 1 elements"))

names(df)
newnames = ["Name","length","rank","Aut_B","Aut_R","Aut_F","Aut_C","Aut_B_stored","Aut_R_stored","Aut_F_stored","Aut_C_stored"]
select!(df,newnames)


=#

