using QuantumAutomorphismGroups
using Oscar
using CSV
using DataFrames
using Statistics
using PlotThemes
using Plots
using StatsPlots


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


function get_unified_df()
    df = read_df(0)

    mat_names = String.(df[!,:Name])

    unified_df = DataFrame(Name=mat_names, deg_bound4=zeros(Int,length(mat_names)), deg_bound10=zeros(Int,length(mat_names)), deg_bound15=zeros(Int,length(mat_names)), no_deg_bound=zeros(Int,length(mat_names)))
    
    unified_df[!,:no_deg_bound] = read_df(-1)[!,:no_deg_bound]
    unified_df[!,:deg_bound4] = read_df(4)[!,:deg_bound4]
    unified_df[!,:deg_bound10] = read_df(10)[!,:deg_bound10]
    unified_df[!,:deg_bound15] = read_df(15)[!,:deg_bound15]
    unified_df[!,:max_deg] = read_df(-1)[!,:max_deg]
\
    datasets =[:no_deg_bound, :deg_bound4, :deg_bound10, :deg_bound15]

    #If the value is -1, we replace it by 1 week in seconds (7*24*60*60)
    for dataset in datasets
        for i in 1:length(unified_df[!,dataset])
            if unified_df[i,dataset] == -1
                unified_df[i,dataset] = 7*24*60*60
            end
        end
    end

    
    return unified_df
end

function print_speed_comparison()
    df = get_unified_df()

    theme(:wong)
    default(fmt=:png,
    size=(800,600),
    legendfontsize=8,
    guidefontsize=8,
    tickfontsize=8,
    titlefontsize=10,
    dpi=300,
    markershape=:cross,
    markercolor=:black,
    markersize=3,
    legend=:topright)

    plot(
    qqplot(df[!,:no_deg_bound],df[!,:deg_bound4],xlabel="nc-buchberger vs degree bound 4",scale=:log10,ylabel="Time (s)"),
    qqplot(df[!,:no_deg_bound],df[!,:deg_bound10],xlabel="nc-buchberger vs degree bound 10",scale=:log10,ylabel="Time (s)"),
    qqplot(df[!,:no_deg_bound],df[!,:deg_bound15],xlabel="nc-buchberger vs degree bound 15",scale=:log10,ylabel="Time (s)"),
    xlims=(0.0001,10^6),
    ylims=(0.0001,10^6))
end

