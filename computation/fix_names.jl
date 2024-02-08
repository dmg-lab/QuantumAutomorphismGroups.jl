# fix names


data_dir = "../data/"
infoFile = ".info"
data_table = "data_table_new.csv"
data_table_old = "data_table_old.csv"

#=
for folder in readdir(data_dir)
    isdir(data_dir * folder) || continue

    for file in readdir(data_dir*folder)
        isfile(data_dir * folder * "/" * file) || continue

        oldname = String(split(file, ".")[1])
        oldM = altNameToMatroid(oldname)
        newname = getName(oldM)
        newM = nameToMatroid(newname)

        bases(oldM) == bases(newM) || error("Bases are not the same for $oldname and $newname")
        
        if oldname != newname
            println("Renaming $oldname to $newname")
            mv(data_dir * folder * "/" * file, data_dir * folder * "/" * newname * infoFile)
        end
    end
end

open(data_dir*data_table,"w") do f
    for line in eachline(data_dir*data_table_old)
        line[1] == 'N' && continue
        oldname = String(split(line, ",")[1])
        oldM = altNameToMatroid(oldname)
        newname = getName(oldM)
        newM = nameToMatroid(newname)

        bases(oldM) == bases(newM) || error("Bases are not the same for $oldname and $newname")
        
        println("Renaming $oldname to $newname")
        println(f, "$newname,$(split(line, ",";limit=2)[2])")

    end
end
=#



