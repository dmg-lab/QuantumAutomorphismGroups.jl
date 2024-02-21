using DataFrames
using ProgressMeter
using CSV


export addAutToDF,
    addToDf,
    loadInfo,
    loadAll,
    loadDf,
    toDataFrame
#=
M = uniform_matroid(0,1)
structure = "bases"
df = addAutToDF(M,true)
select(df,[:Name,:length,:rank, :Aut_B, :LpAut_B])

names(df)
=#
function addAutToDF(M::Matroid, recompute::Bool=false)
    #Load existing DataFrame
    if isfile("../data/data_table.csv")
        df = CSV.read("../data/data_table.csv",DataFrame)
    else
        df = DataFrame(Name=String[])
    end
    df = addAutToDF(M,recompute,df)
    CSV.write("../data/data_table.csv",df)
    return df
end




function addAutToDF(M::Matroid, recompute::Bool=false,df::DataFrame=DataFrame())
    dct = loadDict(M)
    name = getName(M)
    #Load existing DataFrame
    if isfile("../data/data_table.csv")
        df = CSV.read("../data/data_table.csv",DataFrame)
    else
        df = DataFrame(Name=String[])
    end

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

    for key in keys(dct)
        if occursin("Aut",key)
            keywords = split(key,"_") 
            structure = keywords[2]

            #Save initial data
            data["Aut_$(uppercase(structure[1]))_stored"] = true
            if key == "Aut_$(structure)"
                newDataName = "Aut_$(uppercase(structure[1]))"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = isCommutative(dct[key])[1]
                end
            elseif key == "Aut_$(structure)_timed"
                newDataName = "Aut_$(uppercase(structure[1]))_timed"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = dct[key]
                end
            elseif key == "Aut_$(structure)_lp"
                newDataName = "LpAut_$(uppercase(structure[1]))_max_degree"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = maximum(total_degree.(dct[key]))
                end
                maxdeg = newDataName
                newDataName = "LpAut_$(uppercase(structure[1]))"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = isCommutative_lp(dct[key],max(data[maxdeg],2))[1]
                end
            elseif key == "Aut_$(structure)_lp_degree_bound"
                newDataName = "LpAut_$(uppercase(structure[1]))_degree_bound"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = dct[key]
                end
            elseif key == "Aut_$(structure)_lp_timed"
                newDataName = "LpAut_$(uppercase(structure[1]))_timed"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = dct[key]
                end
            end
        end

    end

    addToDf!(data,df)


    return df

end

function addToDf!(data::Dict,df::DataFrame)
    name = data["Name"]
    row = filter(:Name => ==(name),df)
    
    #Only appand if this name does not exist yet
    if size(row)[1] == 0
        append!(df, data, promote=true, cols=:union)
    end
    # Add new columns if they do not exist yet
    new_df = DataFrame(data)
    if Symbol.(names(df)) !== Symbol.(names(new_df))

        for new_col in Symbol.(names(new_df))
            if !(new_col in Symbol.(names(df)))
                println("Adding column $new_col")
                println(typeof(new_df[1,new_col])) 
                insertcols!(df,length(names(df)),new_col => missings(typeof(new_df[1,new_col]), nrow(df)))
            end
        end
        # Add the old columns if they are not in the new dataframe
        for old_col in Symbol.(names(df))
            if !(old_col in Symbol.(names(new_df)))
                println("Adding column $old_col")
                insertcols!(new_df,length(names(new_df)),old_col => missings(typeof(df[1,old_col]), nrow(new_df)))
            end
        end


    end
    new_line = select!(new_df,Symbol.(names(df)))

    df[df.:Name .== name, :] = new_line
    return df
end

#=
d
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

function loadDf()
    CSV.read("../data/data_table.csv",DataFrame)
end

 

#=
info = loadInfo(non_fano_matroid())
split.(keys(info),"_")
I = info["Aut_bases_lp"]

maximum(total_degree.(I)) # To get the maximum degree of the polynomials
AbstractAlgebra.groebner_basis(I)

=#

function loadInfo(name::String,show::Bool=false)

    data_dir = "../data/"
    infoFile = ".info"
    table_dir = data_dir*"data_table.csv"

    folder = split(name,"_")[1]
    
    path = folder * "/"* name
    fullpath = data_dir * path * infoFile

    #Check if folder exists, if not create it
    isdir(data_dir * folder) || throw(ArgumentError("No data for this matroid exists"))


    #Check if info exists, if not create it
    if isfile(fullpath) 
        
        df = CSV.read(table_dir,DataFrame)

        show && show(filter(:Name=> ==(name),df))

        return  loadDict(fullpath)
    else
        throw(ArgumentError("No data for this matroid exists"))
    end

end


loadInfo(M::Matroid) = loadInfo(getName(M))





#=
loadAll(true)
=#
function loadAll(prettyprint::Bool=false)
    data_dir = "../data/"

    data = Dict{String,Dict{String,Any}}()
    
    for folder in readdir(data_dir)
        if isdir(data_dir * folder)
            for file in readdir(data_dir * folder)
                if isfile(data_dir * folder * "/" * file)
                    if prettyprint
                        path = data_dir * folder * "/" * file
                        cmd = pipeline(`jq . $path `,` sponge $path`)
                        run(cmd) 
                    end
                    #cut of the .info
                    fileName = split(file,".")[1]
                    data[fileName] = loadDict(data_dir * folder * "/" * file)
                end
            end
        end
    end
    return data
end


#=
db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]
cursor=Polymake.Polydb.find(collection, Dict("N_ELEMENTS"=>n))
append!(Droids,Matroid.(cursor))

M

computed = []
for (root, dirs, files) in walkdir("../data")
    for file in files
        name, ext = splitext(file)
        ext != ".info" && continue
        println("Processing $name")
        fullpath = joinpath(root, file)
        M = nameToMatroid(name)
        dct = loadDict(M)
        if haskey(dct,"Aut_bases_lp") && r
            push!(computed,name)
        end
    end
end

comp_37 = filter(x->occursin("r3n7",x),computed)


db = Polymake.Polydb.get_db()
collection = db["Matroids.Small"]
cursor=Polymake.Polydb.find(collection, Dict("N_ELEMENTS"=>7,"RANK"=>3))
x = Matroid.(cursor)
computed 


for (root, dirs, files) in walkdir("../data")
    for file in files
        name, ext = splitext(file)
        ext != ".info" && continue
        println("Processing $name")
        fullpath = joinpath(root, file)
        M = nameToMatroid(name)
        dct = loadDict(M)
        if haskey(dct,"Aut_bases_lp")
            dct["Aut_bases_lp_degree_bound"] = length(M)^2 +2
        end
        saveDict(fullpath,dct)
    end
end

M=non_fano_matroid()
length(M)^2+2
loadDict(nameToMatroid(getName(M)))
=#

#= Save and Load File
using Oscar
using DataFrames
using ProgressBars
using CSV


path = "../data/data_table.csv"
df = CSV.read(path,DataFrame)
select!(df,:Aut_B_timed)
the bear
filter(:Name=> ==("r1n7_1"),df)

df = select(df,[:Name,:length,:rank, :Aut_B,:Aut_C,:Aut_B_stored,:Aut_C_stored])
sort!(df,[:length,:rank])
CSV.write(path,df)
loadInfo("r2n6_0003")

=#

#=
path = "../data/data_table.csv"
df = CSV.read(path,DataFrame)
names(df)
=#

#=
i=0
for (root, dirs, files) in walkdir("../data")
    for file in files
        i+=1
    end
end

p = Progress(i)
if isfile("../data/data_table.csv")
    df = CSV.read("../data/data_table.csv",DataFrame)
else
    df = DataFrame(Name=String[])
end

for (root, dirs, files) in walkdir("../data")
    for file in files
        name, ext = splitext(file)
        if ext != ".info" 
            next!(p) 
            continue
        end
        println("Processing $name")
        fullpath = joinpath(root, file)
        M = nameToMatroid(name)
        addAutToDF(M,false,df)
        next!(p)
    end
end

CSV.write("../data/data_table.csv",df)


M= non_fano_matroid()
println.(loadDict(M)["Aut_bases_lp"])
isCommutative_lp(loadDict(M)["Aut_bases_lp"],3)


=#



