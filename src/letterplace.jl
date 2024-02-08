

FreeAssAlgebra = AbstractAlgebra.Generic.FreeAssAlgebra

export lp_groebner_basis, to_FreeAssAlgebra, to_LPRing

@doc raw"""
    to_LPRing(A::FreeAssAlgebra,deg::Int)

Transforms a FreeAssAlgebra to a Singular LPRing with a degree bound `deg`.

# Example
```jldoctest
free, _ = free_associative_algebra(QQ, ["x", "y", "z"])
lp_free, _ = to_LPRing(free, 4)

isa(lp_free, Singular.LPRing)

# output

true
```
"""
to_LPRing(A::FreeAssAlgebra,deg::Int) = Singular.FreeAlgebra(base_ring(A), String.(symbols(A)), deg)

@doc raw"""
    to_FreeAssAlgebra(A::FreeAssAlgebra,a::NCRingElem)

Transforms a Singular NCRingElem to a FreeAssAlgebra element.

# Example
```jldoctest
R, (x, y, z) = Singular.FreeAlgebra(QQ, ["x", "y","z"],6)
free, _ = free_associative_algebra(QQ, ["x", "y", "z"])
f1 = x*y + y*z

F1 = to_FreeAssAlgebra(free, f1)
isa(F1,FreeAssAlgElem)

# output

true
```
"""
function to_FreeAssAlgebra(A::FreeAssAlgebra,a::NCRingElem)
    B = 0
    for (c,e) in zip(Oscar.coefficients(a), Singular.exponent_words(a))
        x = base_ring(A)(c)
        if e == [] 
            B+= x
            continue
        end
        B+= x * (prod(FreeAssAlgElem[A[i] for i in e]))
    end 
    return B
end

@doc raw"""
    lp_groebner_basis(I::Vector{<:FreeAssAlgElem{T}}, n::Int=4; prot::Bool=false) where T<:FieldElem

Computes the Groebner basis of a given ideal I using Letterplace. It returns the Groebner basis and the protocol if `prot` is true


# Example
```jldoctest
M = uniform_matroid(1,2)
rels, u, A = getMatroidRelations(M,:bases)

gb, prot, elapsed = lp_groebner_basis(rels, 4, prot=true)
length(gb)

# output 

4
```
"""
function lp_groebner_basis(I::Vector{<:FreeAssAlgElem{T}},
    n::Int=4;
    prot::Bool=false) where T<:FieldElem 
    length(I) == 0 && return I

    AA_ring = parent(I[1])


    LP_Ring, _   = to_LPRing(AA_ring, n)

    LP_I_gens = LP_Ring.(I)

    I = Singular.Ideal(LP_Ring, LP_I_gens)

    gb = nothing
   
    if !prot 
        gb, elapsed = @timed gens(Singular.std(I)) #Compute the standard basis of an ideal 
        return gb, elapsed
    end

    # If Protocal is requested
    io = IOBuffer()

    old_stdout = stdout
    rd, wr = redirect_stdout()

    Singular.with_prot(true) do; 
        gb, elapsed = @timed gens(Singular.std(I)) #Compute the standard basis of an ideal 
    end

    redirect_stdout(old_stdout)
    close(wr)
    write(io, read(rd))
    close(rd)

    prot = String(take!(io))


    return to_FreeAssAlgebra.(Ref(AA_ring),gb), prot, elapsed
end

