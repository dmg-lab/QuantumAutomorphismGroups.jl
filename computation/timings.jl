using Oscar: FreeAssAlgIdeal
using QuantumAutomorphismGroups
using Oscar
using CSV
using DataFrames
using UnicodePlots
using Statistics

df = CSV.read("../data/data_table.csv",DataFrame)

n = String(df[198,:Name])




cs = df[!,:Aut_C]
cs = [x for x in cs if !ismissing(x)] 

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

