using QuantumAutomorphismGroups
using Oscar
using CSV
using DataFrames
using Statistics
using Distributed



function write_df(df::DataFrame, runner_id::Int)
  CSV.write("../computation/Timing/timings_df_$(runner_id).csv",df)
end

function read_df(runner_id::Int)
  if isfile("../computation/Timing/timings_df_$(runner_id).csv") 
    return CSV.read("../computation/Timing/timings_df_$(runner_id).csv",DataFrame)
  else
    return read_df(0)
  end
end




function compute_all_timings(runner_id::Int=-1, recompute::Bool=false)
  df = read_df(runner_id)
  mat_names = String.(df[!,:Name])
  for name in mat_names
    compute_timing(name,runner_id,runner_id,recompute)
      println("$(name) done, for deg_bound = $(runner_id)")
  end
end



function compute_timing(mat_name::String, deg_bound::Int=-1, runner_id::Int=0, recompute::Bool=false)
  @assert deg_bound in [-1,4,10,15]

  if isfile("../computation/Timing/timings_df_$(runner_id).csv") 
    df = read_df(runner_id)
  else
    df = read_df(0)
  end


  row = df[findfirst(x->x==mat_name,df[!,:Name]),:]
  if deg_bound == 4 
    col = :deg_bound4
  elseif deg_bound == 10
    col = :deg_bound10
  elseif deg_bound == 15
    col = :deg_bound15
  else
    col = :no_deg_bound
  end

  !recompute && row[col] > -1 && return row[col]
  
  
  M = nameToMatroid(mat_name)
  qAut, _ , _  = getMatroidRelations(M,:bases)
  groebner_basis(qAut,deg_bound)
  
  if deg_bound == -1
    gb, elapsed = @timed groebner_basis(qAut)
    deg = reduce(max, map(x -> reduce(max, length.(x.exps)), gb))
    row[:max_deg] = deg
  else
    _, elapsed = @timed groebner_basis(qAut,deg_bound)
  end
  
  row[col] = elapsed

  write_df(df,runner_id)

  return elapsed
end  



#=

compute_all_timings(-1,true)

=#



#=
timings  = df[!,:Aut_B_timed]
# number of missing values
sum(ismissing.(timings))
# number of non missing values
length(timings[.!ismissing.(timings)])
nonmtimings = Float64[ x/60 for x in timings if !ismissing(x)]
#median  in seconds rounded to 2 decimals
round(median(nonmtimings),digits=2)
#mean in minutes rounded to 2 decimals
round(mean(nonmtimings),digits=2)
#std in minutes rounded to 2 decimals
round(std(nonmtimings),digits=2)
#Highest value in hours rounded to 2 decimals
round(maximum(nonmtimings)/3600,digits=2)
#Plot histogram
histogram(nonmtimings,nbins=15,closed=:left)
=#
