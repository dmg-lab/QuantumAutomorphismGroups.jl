using Oscar 

using CSV
using DataFrames
using DataFrames.PrettyTables

path_to_csv = try ARGS[1] catch _  "../computation/dummy.csv" end
path_to_tex = try ARGS[2] catch _ "../computation/test.tex" end

path_to_csv="../data/data_table.csv"
path_to_tex="../../matroid_quantum_automorphism.tex"

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

write(path_to_tex,str)


