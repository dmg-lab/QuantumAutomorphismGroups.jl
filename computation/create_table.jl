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
    
    #hex = map(x->"\texttt{"*String(split(String(x),"_")[2])*"}",df[!,:Name])
    hex = map(x->"\texttt{"*altName(nameToMatroid(String(x)))*"}",df[!,:Name])
    insertcols!(df,2, :hex => hex) 

    #sorting helper
    sh = map(x->1/parse(Int,split(String(x),"_")[2], base=16),df[!,:Name])    
    insertcols!(df,2, :sorting_helper => sh) 
    # All the extra information

    gs = map(x->girth(nameToMatroid(String(x))),df[!,:Name])
    insertcols!(df,:girth => gs)
    ss = map(x->Oscar.order(automorphism_group(nameToMatroid(String(x)))),df[!,:Name])
    insertcols!(df, Symbol("ord(Aut)") => ss)

    dd = map(x-> custom_deg(loadInfo(String(x))["Aut_bases"]),df[!,:Name])
    insertcols!(df, :deg => dd)

    cf = map(x->length(nonbases(nameToMatroid(String(x)))),df[!,:Name])
    insertcols!(df, :n_nonbases => cf)
    sort!(df, [:length,:rank,:sorting_helper])
    select!(df, Not([:sorting_helper]))
    return df
end


s = getName(uniform_matroid(2,3))
altName(uniform_matroid(2,3))


function getHeader(df::DataFrame)

    headerDict = Dict(
    :Name => "Name",
    :hex => "Matroid",
    :Aut_B => latex_cell"$\QAut{\pB}{M}$",
    :Aut_C => latex_cell"$\QAut{\pC}{M}$",
    :n_nonbases => latex_cell"\#nonbases",
    Symbol("ord(Aut)") => latex_cell"$|\Aut(\sM)|$",
    :girth => latex_cell"$\girth(\sM)$",
    :length => latex_cell"$|E(\sM)|$",
    :rank => latex_cell"$\rank(\sM)$",
    :n_loops => latex_cell"\#nloops",
    :CAut => latex_cell"$\Aut(M)$",
    :deg => latex_cell"$d(\pB)$",
    )

    return map(x->headerDict[Symbol(x)],names(df))
end




function toLongtable(df::DataFrame, label::String="Tab:computational-results", caption::String="enter caption here")
    header = getHeader(df)

    io = IOBuffer()
    CSV.write(io,df)
    str2 = String(take!(io))
    str2 = DataFrames.pretty_table(io,df,backend = Val(:latex), header=header)
    str2 = String(take!(io))

    str2 = replace(str2, r"\\hline"m => raw"");
    str2 = replace(str2, r"tabular"m => raw"longtable");

    str2 = replace(str2, r"\\textbf\{(.*?)\}"m => s"\1");
    reg = r"(\\begin\{longtable\}\{r*\}[\r\n]+)([\s^\r\n]+)([^\r\n]+)"m
        replstring = s"""
        \1
          \\caption{cpt10n} \\label{libel} \\\\
          \\toprule
          \3 \\midrule
          \\endfirsthead
          \\caption{cpt10n (continued)} \\\\  
          \3 \\midrule
          \\endhead
          \\endfoot
          \\endlastfoot"""
    str2 = replace(str2, reg => replstring);
str2 = replace(str2, r"\\endhead\n\\endfoot"m => "\\endhead\n\\multicolumn{$(length(header))}{r@{}}{continued on next page \\ldots}\n\\endfoot")
    str2 = replace(str2, r"libel"m => label);
str2 = replace(str2, r"\\end\{longtable\}"m => s"  \\bottomrule\n\\end{longtable}");
    #str2 = replace(str2, r"(?<=\\begin\{longtable\})(\{r*\}\n)(?:.|\n)*?"m=>s"\1 \\hline")
    
    str2 = replace(str2, r"\\begin\{longtable\}\{r*\}\n"m => "\\begin{longtable}{"*"l"*"r"^(length(header)-1)* "}");
     
    #delete last \hline
    str2 = replace(str2, r"\\\\ \\hline\n\\label" => "\n\\label");
    str2 = replace(str2, r"cpt10n"m => caption);

    str2 = replace(str2, r"\\textbackslash\{\}"m => "\\");
    str2 = replace(str2, r"\\{"m => s"{");
    str2 = replace(str2, r"\\}"m => s"}");
     
    #Checkmark means it has a non commutative quantum automorphism group
    str2 = replace(str2, r"false"m => raw"$\times$");
    str2 = replace(str2, r"true"m => raw"$\checkmark$");
    str2 = replace(str2, r"missing"m => raw"?");
    
    return str2
end


function putDataInLateX(path_to_csv,path_to_tex)
    caption = match(r"(?<=\\caption\{)(?:.)*?(?=\\label)",str).match
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



function replaceByLabel(str1::String, label::String, df::DataFrame)
        longtab(label) = Regex(s"(?:\\begin\{longtable\}\{[^\r\n]*\}\n(?:.)*\\label\{)"*label*s"(?:(?!\\begin\{longtable\}).|\n)*?\\end\{longtable\}")
        re = longtab(label)
        m = match(re,str1)
        isnothing(m) && return str1
        oldtable = m.match
        m2 = match(r"(?<=\\caption\{)(?:.)*?(?=\} \\label)",oldtable)
        isnothing(m2) && return str1
        caption = String(m2.match)
        str2 = toLongtable(df,label,caption)
        str = replace(str1, re=>str2)
        return str 
end
#=
csv="../data/data_table.csv"
tex="../../matroid_quantum_automorphism.tex"

df = getExtraInf(csv)
print(toLongtable(df))
str = read(tex, String);
print(str)

str = replaceByLabel(str,"Tab:computational-results-1",df)
print(str)


csv="../data/data_table.csv"
df = getExtraInf(csv)
=#

function replaceAll()

    csv="../data/data_table.csv"
    tex="../../paper.tex"

    df = getExtraInf(csv)
    str = read(tex, String);
    df_nm=filter(x->!ismissing(x[:Aut_B]) && !ismissing(x[:Aut_C]),df)
    #Table in with both Aut_B and Aut_C are false
    df1 = filter(x->x[:Aut_B]==false && x[:Aut_C]==false,df_nm)
    df1 = select(df1,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)"),:deg])
    newstr = replaceByLabel(str,"Tab:computational-results-1",df1)


    #Table in with both Aut_B and Aut_C are true
    df2 = filter(x->x[:Aut_B]==true && x[:Aut_C]==true,df_nm)
    df2 = select(df2,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)"),:deg])
    newstr = replaceByLabel(newstr,"Tab:computational-results-2",df2)


    #Table in with Aut_B is true and Aut_C is false
    df3 = filter(x->x[:Aut_B]==true && x[:Aut_C]==false,df_nm)
    df3 = select(df3,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)"),:deg])
    newstr = replaceByLabel(newstr,"Tab:computational-results-3",df3)

    #Table in with Aut_B is true and Aut_C is false
    df4 = filter(x->x[:Aut_B]==false && x[:Aut_C]==true,df_nm)
    df4 = select(df4,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)"),:deg])
    newstr = replaceByLabel(newstr,"Tab:computational-results-4",df4)

    #Table in with Aut_B is true and Aut_C is missing
    df5 = filter(x->x[:Aut_B]==true && ismissing(x[:Aut_C]),df)
    df5 = select(df5,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)"),:deg])
    newstr = replaceByLabel(newstr,"Tab:computational-results-5",df5)

    #Table in with Aut_B is true and Aut_C is missing
    df6 = filter(x->x[:Aut_B]==false && ismissing(x[:Aut_C]),df)
    df6 = select(df6,[:hex,:length,:rank,:girth,:n_nonbases,Symbol("ord(Aut)"),:deg])
    newstr = replaceByLabel(newstr,"Tab:computational-results-6",df6)


    write(tex,newstr);
end

#=
replaceAll()
=#

