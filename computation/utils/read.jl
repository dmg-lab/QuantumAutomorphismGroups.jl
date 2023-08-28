using FileWatching


function readInput(filename::String="input.txt")
    inpt_D = Dict{String, Any}()

    inpt = filter(!isempty,readlines(filename))
    for i in eachindex(inpt)
        if inpt[i][1] == '-' && inpt[i][14] == '0'
            name = inpt[i+1]
            M = inpt[i+2]
            inpt_D[name] = M
        end
    end
    return  inpt_D
end

In_loc =  try ARGS[1] catch _ "input.txt" end
ipt = readInput(In_loc)
for (k,v) in ipt
    println(k*":"*v)
end



