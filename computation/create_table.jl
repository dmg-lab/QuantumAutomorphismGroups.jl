using Oscar 

using QuantumAutomorphismGroups

using CSV
using DataFrames
using DataFrames.PrettyTables

csv = try ARGS[1] catch _  "../computation/dummy.csv" end
tex = try ARGS[2] catch _ "../computation/test.tex" end


function toString(v:: Vector{String})
    return "[" * join(v, ", ") * "]"
end


function toString(v:: Vector{Vector})
    return map(x -> reduce(*, string.(x)), v) |> cs -> filter(x-> x!= "", cs) |> toString
end

function toString(v:: Vector{Vector{Int}})
    return map(x -> reduce(*, string.(x)), v) |> cs -> filter(x-> x!= "", cs) |> toString
end


function getExtraInf(path_to_csv)
    df = CSV.read(csv,DataFrame)
    names(df)
    df = select(df,[:Name,:length,:rank,:Aut_B,:Aut_C])
    
    hex = map(x->split(String(x),"_")[2],df[!,:Name])
    insertcols!(df,2, :hex => hex) 
    # All the extra information

    gs = map(x->girth(nameToMatroid(String(x))),df[!,:Name])
    insertcols!(df,:girth => gs)
    ss = map(x->Oscar.order(automorphism_group(nameToMatroid(String(x)))),df[!,:Name])
    insertcols!(df, Symbol("ord(Aut)") => ss)

    cf = map(x->length(nonbases(nameToMatroid(String(x)))),df[!,:Name])
    insertcols!(df, :n_nonbases => cf)
    sort!(df, [:length,:rank,:n_nonbases,:girth,:hex])
    return df
end

function getHeader(df::DataFrame)

    headerDict = Dict(
    :Name => "Name",
    :hex => "Name",
    :Aut_B => latex_cell"$\QAut{\pB}{M}$",
    :Aut_C => latex_cell"$\QAut{\pC}{M}$",
    :n_nonbases => latex_cell"n\_nonbases",
    Symbol("ord(Aut)") => latex_cell"$\ord(\Aut(M))$",
    :girth => latex_cell"$\girth(M)$",
    :length => latex_cell"$|E(M)|$",
    :rank => latex_cell"$r(M)$",
    :n_loops => latex_cell"n\_loops",
    :CAut => latex_cell"$\Aut(M)$",

    )

    return map(x->headerDict[Symbol(x)],names(df))
end



function toLongtable(df::DataFrame, label::String="Tab:computational-results")
    header = getHeader(df)

    io = IOBuffer()
    CSV.write(io,df)
    str2 = String(take!(io))
    str2 = DataFrames.pretty_table(io,df,backend = Val(:latex), header=header)
    str2 = String(take!(io))

    str2 = replace(str2, r"\\hline"m => raw"");
    str2 = replace(str2, r"tabular"m => raw"longtable");

    str2 = replace(str2, r"\\\\\\\\"m => "\\\\ \\hline");


str2 = replace(str2, r"(?:.)*?(?=\\end\{longtable})"m => "\\label{"*label*"} \\\\ \\hline\n")
    str2 = replace(str2, r"(?<=\\begin\{longtable\})(\{r*\}\n)(?:.|\n)*?"m=>s"\1 \\hline")
    str2 = replace(str2, r"\\begin\{longtable\}\{r*\}\n"m => "\\begin{longtable}{"*"|l|"*"r|"^(length(header)-1)* "}\n");
    
    #delete last \hline
    str2 = replace(str2, r"\\\\ \\hline\n\\label" => "\n\\label");
     
    #Checkmark means it has a non commutative quantum automorphism group
    str2 = replace(str2, r"false"m => raw"$\times$");
    str2 = replace(str2, r"true"m => raw"$\checkmark$");
    str2 = replace(str2, r"missing"m => raw"?");
    
    return str2
end


function putDataInLateX(path_to_csv,path_to_tex)

    df = CSV.read(path_to_csv,DataFrame)
    nams =[:Name,:length,:rank, :Aut_B,:Aut_C]
    df = select(df,nams)



    str2 = toLongtable(df)

    str = read(path_to_tex, String);


    re = r"(?<=\\small\n)((?:.|\n)*?)(?=\\end{center})"

    str = replace(str, re=>"$(str2)\n")

    write(path_to_tex,str);
    println("updated $path_to_tex")

    return 
end

function replaceByLabel(str1::String, label::String, str2::String)
    longtab(label) = Regex(raw"(?:\\begin\{longtable\})(?:(?!\\begin\{longtable\}).|\n)*?\\label\{"*label*raw"\} \\\\ \\hline\n\\end{longtable}")
        re = longtab(label)
        str = replace(str1, re=>str2)
        return str 
end



csv="../data/data_table.csv"
tex="../../matroid_quantum_automorphism.tex"

df = getExtraInf(csv)
str = read(tex, String);

df_nm=filter(x->!ismissing(x[:Aut_B]) && !ismissing(x[:Aut_C]),df)
#Table in with both Aut_B and Aut_C are false
df1 = filter(x->x[:Aut_B]==false && x[:Aut_C]==false,df_nm)
df1 = select(df1,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)")])
newstr = replaceByLabel(str,"Tab:computational-results-1",toLongtable(df1,"Tab:computational-results-1"))


#Table in with both Aut_B and Aut_C are true
df2 = filter(x->x[:Aut_B]==true && x[:Aut_C]==true,df_nm)
df2 = select(df2,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)")])
newstr = replaceByLabel(newstr,"Tab:computational-results-2",toLongtable(df2,"Tab:computational-results-2"))


#Table in with Aut_B is true and Aut_C is false
df3 = filter(x->x[:Aut_B]==true && x[:Aut_C]==false,df_nm)
df3 = select(df3,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)")])
newstr = replaceByLabel(newstr,"Tab:computational-results-3",toLongtable(df3,"Tab:computational-results-3"))

#Table in with Aut_B is true and Aut_C is false
df4 = filter(x->x[:Aut_B]==false && x[:Aut_C]==true,df_nm)
df4 = select(df4,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)")])
newstr = replaceByLabel(newstr,"Tab:computational-results-4",toLongtable(df4,"Tab:computational-results-4"))

#Table in with Aut_B is true and Aut_C is missing
df5 = filter(x->x[:Aut_B]==true && ismissing(x[:Aut_C]),df)
df5 = select(df5,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)")])
newstr = replaceByLabel(newstr,"Tab:computational-results-5",toLongtable(df5,"Tab:computational-results-5"))

write(tex,newstr);

