using DataFrames
using ProgressMeter
using CSV


export addAutToDF,
    addToDf,
    loadInfo,
    loadAll,
    loadDf,
    toDataFrame,
    getVanishingVariables,
    isIncluded,
    getSubgroups,
    subgroupDict,
    getClasses,
    getPath,
    addToInfo

#=
M = uniform_matroid(0,1)
structure = "bases"
df = addAutToDF(M,true)
loadInfo(fano_matroid())
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
    df = addAutToDF(M,df,recompute)
    CSV.write("../data/data_table.csv",df)
    return df
end




function addAutToDF(M::Matroid, df::DataFrame, recompute::Bool=false)
    dct = loadDict(M)
    name = getName(M)
    #Load existing DataFrame

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
                newDataName = "LpAut_$(uppercase(structure[1]))_n_vvariables"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = length(getVanishingVariables(dct[key]))
                end
                maxdeg = data["LpAut_$(uppercase(structure[1]))_max_degree"]
                newDataName = "LpAut_$(uppercase(structure[1]))"
                if recompute || !haskey(data, newDataName) || ismissing(data[newDataName])
                    data[newDataName] = isCommutative_lp(dct[key],max(maxdeg,2))[1]
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
getPath(non_fano_matroid())
=#
function getPath(name::String)
    data_dir = "../data/"
    dir, _ = split(name,"_")
    return data_dir* dir*"/" * name * ".info"
end
getPath(M::Matroid) = getPath(getName(M)) 

#=
info = loadInfo(fano_matroid())
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
        return  load(fullpath)
    else
        throw(ArgumentError("No data for this matroid exists"))
    end

end

loadInfo(M::Matroid) = loadInfo(getName(M))

#=
x, timimg = @timed loadAll(false,false);
=#
function loadAll(resave::Bool=false, prettyprint::Bool=true, showprogress::Bool=true)
    data_dir = "../data/"

    data = Dict{String,Dict{String,Any}}()
    

    filenames = []
    for (root, _, files) in walkdir(data_dir)
        for file in files
            name, ext = splitext(file)
            if ext != ".info" 
                continue
            end
            path = root * "/" * file
            push!(filenames,(name,path))
        end
    end
    @info "$(length(filenames)) files to process."

    if showprogress
        @showprogress for (name,path) in filenames 
            dct = load(path)
            data[name] = dct
            if resave
                save(path,dct)
            end
            if prettyprint
                cmd = pipeline(`jq . $path `,`sponge $path`)
                run(cmd)
            end
        end
    else
        for (name,path) in filenames 
            dct = load(path)
            data[name] = dct
            if resave
                save(path,dct)
            end
            if prettyprint
                cmd = pipeline(`jq . $path `,`sponge $path`)
                run(cmd)
            end
        end
    end

    return data
end

#=
M = non_fano_matroid()
=#
function addToInfo(name::String, key::String, value::Any)
    path = getPath(name)
    @assert isfile(path) "No data for this matroid exists"
    info = loadInfo(name)
    info[key] = value
    save(path,info)
    return info
end
addToInfo(M::Matroid, key::String, value::Any) = addToInfo(getName(M),key,value)

#=
M = non_fano_matroid()
info = loadInfo(M)

x = length(getVanishingVariables(M))
=#

function getVanishingVariables(M::Matroid) 
    info = loadInfo(M)
    haskey(info,"Aut_bases_lp") || return nothing
    I = info["Aut_bases_lp"]
    return getVanishingVariables(I)
end

function getVanishingVariables(I::Vector{<:FreeAssAlgElem{T}}) where T <:FieldElem
    deg = maximum(total_degree.(I))
    Alg = parent(I[1])
    
    LP_Ring, _  = to_LPRing(Alg, deg)
    I = Singular.Ideal(LP_Ring, LP_Ring.(I))
    
    ans = []
    for var in gens(Alg)
        iszero(reduce(LP_Ring(var),I)) && push!(ans,var)
    end
    return ans
end


#=
M = uniform_matroid(3,4)
isIncluded(M,M)
=#

# Check weather the Quantum automorphism group of M is included in the Quantum automorphism group of N
function isIncluded(M::Matroid, N::Matroid)
    @assert length(M) == length(N)
    infoM = loadInfo(M)
    infoN = loadInfo(N)
    haskey(infoM,"Aut_bases_lp") || return nothing
    haskey(infoN,"Aut_bases_lp") || return nothing
    I = infoM["Aut_bases_lp"]
    J = infoN["Aut_bases_lp"]
    return isIncluded(J,I) #Switched the order
end


# Check weather the Ideal I is included in the Ideal of J
function isIncluded(I::Vector{<:FreeAssAlgElem{T}}, J::Vector{<:FreeAssAlgElem{T}}) where T <:FieldElem
    deg = max(maximum(total_degree.(I)),maximum(total_degree.(J)))
    Alg = parent(I[1])
    
    LP_Ring, _  = to_LPRing(Alg, deg)
    J = Singular.Ideal(LP_Ring, LP_Ring.(J))
    for func in I
        if !iszero(reduce(LP_Ring(func),J))
            return false
        end
    end 
    return true
end

#=
M = uniform_matroid(2,4)
sg = getSubgroups(M)

N = non_fano_matroid()
getSubgroups(N)


QuantumAutomorphismGroups.addSubgroupsToInfo()
=#
function addSubgroupsToInfo()
    dct = loadAll(false,false,false)
    @showprogress for (name,info) in pairs(dct)
        M = nameToMatroid(name)
        subgroups = getSubgroups(M)
        isnothing(subgroups) && continue
        addToInfo(M,"Aut_bases_lp_subgroups",subgroups)
        
    end
end




#=
using Oscar
using DataFrames
using CSV

path = "../data/data_table.csv"
df = CSV.read(path,DataFrame)
sort!(df,:Aut_B_timed)
mat_names = String.(df[!,:Name])
unique!(mat_names)





df = select(df,[:Name,:length,:rank, :Aut_B,:Aut_C,:Aut_B_stored,:Aut_C_stored])
sort!(df,[:length,:rank])
CSV.write(path,df)
loadInfo("r2n6_16ef")
r3n6_001ff
H = nameToMatroid("r2n6_16ef")
H = matroid_from_bases([[2, 4, 6], [3, 4, 6], [2, 5, 6], [3, 5, 6], [4, 5, 6]],6)
H1 = matroid_from_bases([[2, 3], [2, 4], [3, 4], [2, 5], [3, 5], [4, 5], [2, 6], [3, 6], [4, 6], [5, 6]],6)
bases(H)

=#

#=
df = loadDf()
names(df)
df1 = select(df,[:Name])
df1 = select(df,[:Name,:length,:rank, :LpAut_B, :Aut_B])
print(df1)
sort!(df1,:LpAut_B_max_degree)
sort!(df1,:LpAut_B_timed)
getName(fano_matroid())
getName(non_fano_matroid())
info = loadInfo(fano_matroid())

gb_lp = info["Aut_bases_lp"]
gb = AbstractAlgebra.groebner_basis(gb_lp)

M = nameToMatroid("r3n7_000007fff")
cyclic_flats(non_fano_matroid())
print(df1)
=#
#=
info = loadInfo(non_fano_matroid())
info["Aut_bases_lp"]

getVanishingVariables(non_fano_matroid())

getVanishingVariables(nameToMatroid("r3n7_3f7efffff"))

pop_first!(cyclic_flats(M))
cyclic_flats(nameToMatroid("r3n7_3f7eefd7f"))
getSubgroups(nameToMatroid("r3n7_3f7eefd7f"))
=#
