using Oscar
using Polymake
using Combinatorics


function matroidIncidence(M::Matroid,t::Symbol=:circuits)
    rows = eval(t)(M);

    groundset = matroid_groundset(M);
    incidence = zeros(Int64,length(rows),length(groundset));

    for i in 1:length(rows)
        for j in 1:length(groundset)
            if groundset[j] in rows[i]
                incidence[i,j] = 1;
            end
        end
    end
    
    return incidence;
end


function toAdjacency(Inc::Matrix{Int64})
    t = hcat(zeros(Int64,size(Inc,1),size(Inc,1)),Inc)
    t1 = hcat( transpose(Inc),zeros(Int64,size(Inc,2),size(Inc,2)))
    adj = vcat(t,t1)
    return adj
end


function Base.contains(V::Vector,v::Vector)
    isempty(v)  && return false
    S = string(V)
    s = string(v) 
    x = findfirst("[",s)
    s = s[x[1]+1:end-1]
    return contains(S,s)
end

function deleteContaining!(V::Vector,verbose::Bool=false)
     
    isempty(V) && return V
    for size in 1:length(V)
        i = 1 
        while i <= length(V[size])
            for bigger in size+1:length(V)
                isempty(V[bigger]) && continue
                j = 1
                while j <= length(V[bigger])
                    if  contains(V[bigger][j],V[size][i]) 
                        splice!(V[bigger],j) 
                    else
                        j+=1
                    end
                end
            end
        i+=1
        end
    end
    return V
end
