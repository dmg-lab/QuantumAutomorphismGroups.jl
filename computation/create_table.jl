using Oscar 

using QuantumAutomorphismgGroups

using CSV
using DataFrames
using DataFrames.PrettyTables

function putDataInLateX(path_to_csv,path_to_tex)

    df = CSV.read(path_to_csv,DataFrame)
    nams =[:Name,:length,:rank, :Aut_B,:Aut_C]
    df = select(df,nams)

    nams2 = ["Name","length","rank",
    latex_cell"$\QAut{\pB}{M}$",
    latex_cell"$\QAut{\pC}{M}$"]


    io = IOBuffer()
    CSV.write(io,df)
    str2 = String(take!(io))
    str2 = DataFrames.pretty_table(io,df,backend = Val(:latex), header=nams2)
    str2 = String(take!(io))

    str2 = replace(str2, r"\\hline"m => raw"");
    str2 = replace(str2, r"missing"m => raw"?");
    str2 = replace(str2, r"tabular"m => raw"longtable");

    str2 = replace(str2, r"\\\\\\\\"m => "\\\\ \\hline");


    str2 = replace(str2, r"(?:.)*?(?=\\end\{longtable})"m => s"\\label{Tab:computational-results} \n")
    str2 = replace(str2, r"(?<=\\begin\{longtable\})(\{r*\}\n)(?:.|\n)*?"m=>s"\1 \\hline")


    str = read(path_to_tex, String);


    re = r"(?<=\\small\n)((?:.|\n)*?)(?=\\end{center})"

    str = replace(str, re=>"$(str2)\n")

    write(path_to_tex,str);
    println("updated $path_to_tex")

    return 
end

csv = try ARGS[1] catch _  "../computation/dummy.csv" end
tex = try ARGS[2] catch _ "../computation/test.tex" end


#csv="../data/data_table.csv"
#tex="../../matroid_quantum_automorphism.tex"

putDataInLateX(csv,tex)


df = CSV.read(csv,DataFrame)

df = select(df,[:Name,:length,:rank])
df[1,:Name]
M = nameToMatroid(String(df[1,:Name]))
length(loops(M))
girth(M)


function toString(v:: Vector{String})
    return "[" * join(v, ", ") * "]"
end
CS(M) = circuits(M) |> cs -> filter(x->length(x)==2, cs) |> cs -> map(x -> "$(x[1])$(x[2])", cs) |> toString

Oscar.describe(automorphism_group(M))


function toString(v:: Vector{Vector})
    return map(x -> reduce(*, string.(x)), v) |> cs -> filter(x-> x!= "", cs) |> toString
end

function toString(v:: Vector{Vector{Int}})
    return map(x -> reduce(*, string.(x)), v) |> cs -> filter(x-> x!= "", cs) |> toString
end

toString(cyclic_flats(M))
cf = map(x->toString(cyclic_flats(nameToMatroid(String(x)))),df[!,:Name])
insertcols!(df, 4, :cyclic_flats => cf)
cf = map(x->CS(nameToMatroid(String(x))),df[!,:Name])
insertcols!(df, 4, :Symbol("2-circuits") => cf)

gs = map(x->girth(nameToMatroid(String(x))),df[!,:Name])
insertcols!(df, 4, :girth => gs)
ls = map(x->length(loops(nameToMatroid(String(x)))),df[!,:Name])
insertcols!(df, 4, :n_loops => ls)

