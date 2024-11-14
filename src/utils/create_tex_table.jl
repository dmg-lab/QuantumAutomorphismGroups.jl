using DataFrames.PrettyTables

export _get_extra_info, _to_longtable_content


function toString(v:: Vector{String})
    return "[" * join(v, ", ") * "]"
end


function toString(v:: Vector{Vector})
    return map(x -> reduce(*, string.(x)), v) |> cs -> filter(x-> x!= "", cs) |> toString
end

function toString(v:: Vector{Vector{Int}})
    return map(x -> reduce(*, string.(x)), v) |> cs -> filter(x-> x!= "", cs) |> toString
end


function _get_extra_info(path_to_csv="../data/data_table.csv")
    df = CSV.read(path_to_csv,DataFrame)
    names(df)
    df = select(df,["Name","Aut_bases_-1","Aut_circuits_-1","Aut_bases_-1_max_deg"])
    
    hex = map(x->"\\texttt{"*String(split(String(x),"_")[2])*"}",df[!,:Name])
    insertcols!(df,2, :hex => hex) 

    #sorting helper
    sh = map(x->1/parse(Int,split(String(x),"_")[2], base=16),df[!,:Name])    
    insertcols!(df,2, :sorting_helper => sh) 
    # All the extra information

    rk = map(x->rank(matroid_from_matroid_hex(String(x))),df[!,:Name])
    insertcols!(df,:rank => rk)
    
    lg = map(x->length(matroid_from_matroid_hex(String(x))),df[!,:Name])
    insertcols!(df,:length => lg)


    gs = map(x->girth(matroid_from_matroid_hex(String(x))),df[!,:Name])
    insertcols!(df,:girth => gs)
    ss = map(x->Oscar.order(automorphism_group(matroid_from_matroid_hex(String(x)))),df[!,:Name])
    insertcols!(df, Symbol("ord(Aut)") => ss)

    #dd = map(x-> custom_deg(loadInfo(String(x))["Aut_bases"]),df[!,:Name])
    #insertcols!(df, :deg => dd)

    cf = map(x->length(nonbases(matroid_from_matroid_hex(String(x)))),df[!,:Name])
    insertcols!(df, :n_nonbases => cf)
    sort!(df, [:length,:rank,:sorting_helper])
    select!(df, Not([:sorting_helper]))
    
    return df

end



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

function _to_longtable_content(path_to_csv="../data/data_table.csv")
    df = _get_extra_info(path_to_csv)

    df = DataFrame(filter(row -> (!ismissing(row["Aut_bases_-1"]) || !ismissing(row["Aut_circuits_-1"])), eachrow(df)))
    


    ncb_and_ncc = filter(row -> !ismissing(row["Aut_bases_-1"]) && row["Aut_bases_-1"] == false && !ismissing(row["Aut_circuits_-1"]) && row["Aut_circuits_-1"] == false, eachrow(df))
    ncb_and_cc = filter(row -> !ismissing(row["Aut_bases_-1"]) && row["Aut_bases_-1"] == false && !ismissing(row["Aut_circuits_-1"]) && row["Aut_circuits_-1"] == true, eachrow(df))
    ncb_and_mc = filter(row -> !ismissing(row["Aut_bases_-1"]) && row["Aut_bases_-1"] == false && ismissing(row["Aut_circuits_-1"]), eachrow(df))

    #cb_and_ncc is empty?!
    cb_and_cc = filter(row -> !ismissing(row["Aut_bases_-1"]) && row["Aut_bases_-1"] == true && !ismissing(row["Aut_circuits_-1"]) && row["Aut_circuits_-1"] == true, eachrow(df))
    cb_and_ncc = filter(row -> !ismissing(row["Aut_bases_-1"]) && row["Aut_bases_-1"] == true && !ismissing(row["Aut_circuits_-1"]) && row["Aut_circuits_-1"] == false, eachrow(df))
    cb_and_mc = filter(row -> !ismissing(row["Aut_bases_-1"]) && row["Aut_bases_-1"] == true && ismissing(row["Aut_circuits_-1"]), eachrow(df))

    mb_and_cc = filter(row -> ismissing(row["Aut_bases_-1"]) && !ismissing(row["Aut_circuits_-1"]) && row["Aut_circuits_-1"] == true, eachrow(df))
    #does not exist
    mb_and_ncc = filter(row -> ismissing(row["Aut_bases_-1"]) && !ismissing(row["Aut_circuits_-1"]) && row["Aut_circuits_-1"] == false, eachrow(df))
    #select!(df, ["hex", "length", "rank", "girth", "n_nonbases", "ord(Aut)", "Aut_bases_-1_max_deg"])

    list_of_tables = DataFrame.([ncb_and_ncc, ncb_and_cc, ncb_and_mc, cb_and_cc, cb_and_ncc, cb_and_mc, mb_and_cc, mb_and_ncc])
    


    
    #Write all of those to csv files
    for (table, name) in zip(list_of_tables, ["ncb_and_ncc", "ncb_and_cc", "ncb_and_mc", "cb_and_cc", "cb_and_ncc", "cb_and_mc", "mb_and_cc", "mb_and_ncc"])
        table =  select(table, ["hex", "length", "rank", "girth", "n_nonbases", "ord(Aut)", "Aut_bases_-1_max_deg"])
        CSV.write("../notebooks/tables/" * name * ".csv", table)
        #Delete first row of that file
        #run(`sed -i '1d' ../notebooks/tables/$name.csv`)
        #delete everything before the first comma including the comma
        #run(`sed -i 's/^[^,]*,//' ../notebooks/tables/$name.csv`)
        #replace each "," with " & "
        run(`sed -i 's/,/ \& /g' ../notebooks/tables/$name.csv`)
        #replace each newline with " \\" and newline
        run(`sed -i 's/$/ \\\\/' ../notebooks/tables/$name.csv`)
    end


    return df
end


#=
rvlx = min_revlex_basis_encoding(uniform_matroid(1,4))
 r,n = rank(M), length(M) 


  v = Oscar._revlex_basis_to_vector(rvlx)

v1 = zeros(Int, 4 - length(v)%4)
v2 = vcat(v1, v)

=#
