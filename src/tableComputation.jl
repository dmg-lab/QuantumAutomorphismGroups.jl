using DataFrames
using ProgressBars
using CSV


export addAutToDF,
    addToDf,
    loadInfo,
    loadAll,
    loadDf,
    toDataFrame
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
    
    pbar = ProgressBar(total=length(dt)) 
    for (name, dict) in dt
        update(pbar)
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
                data["Aut_$(uppercase(structure[1]))_stored"] = true
                newDataName = "Aut_$(uppercase(structure[1]))"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    println("Checking if $name for $structure is commutative by algebraic means")
                    data[newDataName] = isCommutative(dict[key])[1]
                end
                #newDataName = "Aut_$(uppercase(structure[1]))_extra"
                #if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                #    println("Checking if $name for $structure is commutative using extra techniques")
                #    data[newDataName] = isCommutativeExtra(dict[key])[1]
                #end
                if occursin("timed",key)
                    newDataName = "Aut_$(uppercase(structure[1]))_timed"
                    if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                        println("Adding timing for $structure on $name.")
                        data[newDataName] = dict[key]
                    end
                end

            end 

        end
        #only appand if this name does not exist yet
        if size(row)[1] == 0
            append!(df,data,promote=true,cols=:union)
            continue
        end
        new_df = DataFrame(data)
        if Symbol.(names(df)) !== Symbol.(names(new_df))
            for name in Symbol.(names(new_df))
                if !(name in Symbol.(names(df)))
                    insertcols!(df,length(names(df)),name => missings(typeof(new_df[1,name]),nrow(df)))
                end
            end
        end
        new_line = select!(new_df,Symbol.(names(df)))

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

function loadDf()
    CSV.read("../data/data_table.csv",DataFrame)
end

 

#=
loadInfo(uniform_matroid(1,3))

I = loadInfo("r1n2_1")
I["Aut_bases"]
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




#= Save and Load File
using Oscar
using DataFrames
using ProgressBars
using CSV


path = "../data/data_table.csv"
df = CSV.read(path,DataFrame)
sort!(df,:Aut_B_timed)
select!(df,[:Name,:length,:rank, :Aut_B,:Aut_C,:Aut_B_timed])
println(df)


filter(:Name=> ==("r1n7_1"),df)


df = select(df,[:Name,:length,:rank, :Aut_B,:Aut_C,:Aut_B_stored,:Aut_C_stored])
sort!(df,[:length,:rank])
CSV.write(path,df)
loadInfo("r2n6_16ef")
r3n6_001ff
H = nameToMatroid("r2n6_16ef")
H = matroid_from_bases([[2, 4, 6], [3, 4, 6], [2, 5, 6], [3, 5, 6], [4, 5, 6]],6)
H1 = matroid_from_bases([[2, 3], [2, 4], [3, 4], [2, 5], [3, 5], [4, 5], [2, 6], [3, 6], [4, 6], [5, 6]],6)
bases(H)

print("[")
for b in bases(H)
print("$(b), ")
end
println("]")
@time 
=#

#=
path = "../data/data_table.csv"
df = CSV.read(path,DataFrame)
names(df)
=#

#=

=#



