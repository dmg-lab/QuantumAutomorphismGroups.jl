
using Oscar


FreeAssAlgebra = AbstractAlgebra.Generic.FreeAssAlgebra

# freeAlgebra and makeLetterplaceRing are the same?

#=
strng = """
LIB "freegb.lib";
ring r = 0,(x,y,z),dp;
def R = freeAlgebra(r, 4);  // degree (length) bound 4; the ordering will be degree right lex
setring R;
option(redSB);
option(redTail);
option(prot);  // let us activate the protocol
ideal I = x*y + y*z, x*x + x*y - z; // a non-graded ideal
ideal J = letplaceGBasis(I);
"""

Singular.call_interpreter(strng);

ring, dct = Singular.lookup_library_symbol("Top", "R")
SI = dct[:J]

dct[:J]
=#





function groebner_basis_new(I, prot::Bool=false)
    gb = nothing
   
    if !prot 
        gb = gens(Singular.std(I)) #Compute the standard basis of an ideal 
        return gb
    end

    # If Protocal is requested
    io = IOBuffer()


    old_stdout = stdout
    rd, wr = redirect_stdout()

    Singular.with_prot(true) do; 
        gb = gens(Singular.std(I)) #Compute the standard basis of an ideal 
    end

    redirect_stdout(old_stdout)
    close(wr)
    write(io, read(rd))
    close(rd)

    prot = String(take!(io))


    return prot,gb
end

#=
R, (x, y, z) = Singular.FreeAlgebra(QQ, ["x", "y","z"],6)

f1 = x*y + y*z
f2 = x*x + x*y - z

I = Singular.Ideal(R, [f1, f2])

prot, gb = groebner_basis_new(I)
print(prot)
=#


# Transform FreeAssAlgebra to LetterplaceRing

FreeAlgebra(A::FreeAssAlgebra,deg::Int) = Singular.FreeAlgebra(base_ring(A), String.(symbols(A)), deg)


M = uniform_matroid(1,2)
rels, u, A = getMatroidRelations(M,:bases)
sA, gns   = FreeAlgebra(A, 4)
srels = sA.(rels)

I = Singular.Ideal(sA, srels)

prot, gb = groebner_basis_new(I, true)
gb
print(prot)
# isa(sA,Singular.LPRing)


y = gb[1]
collect(exponents(y))
collect(coefficients(y))


isa(gb[1],NCRingElem)
base_ring(A)
function (A::FreeAssAlgebra)(a::NCRingElem)
    B = 0
    println(a)
    for (c,e) in zip(collect(coefficients(a)), collect(exponents(a)))
        println(FreeAssAlgElem[A[i]^j for (i,j) in enumerate(e) if j != 0])
        B+= base_ring(A)(c)*(prod(FreeAssAlgElem[A[i]^j for (i,j) in enumerate(e)]))
    end 
end

A(gb[1])
