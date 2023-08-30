
Pkg.activate("../")

using Oscar
using Oscar.JSON


include("./utils/save.jl")
include("../src/quantumMatroid.jl")

auts = loadAhoCorasick("./data/r3n7/fano_matroid/Aut_bases.aho")
info = loadDict("./data/r3n7/fano_matroid/stats.info")

gns = info["Aut_bases"][1]
dur = info["Aut_bases"][2]
Alg = parent(gens[1])




    


Alg, _, gns, auts = getIdeal(uniform_matroid(2,3))

isInIdeal(gns,gns,auts)


str = JSON.json(auts)
newAuts = strToAho(str)

saveAho("./tst.tst",auts)
newAuts = loadAho("./tst.tst")

isInIdeal(gns,gns,newAuts)



