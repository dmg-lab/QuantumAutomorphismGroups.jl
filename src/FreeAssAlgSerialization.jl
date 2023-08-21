import Oscar.save_internal
import Oscar.load_internal
import Oscar.encodeType
import Oscar.registerSerializationType
import Oscar.@registerSerializationType
import Oscar.PolyRingElem

using Oscar
import AbstractAlgebra
import AbstractAlgebra.Generic


# Free associative algebra serialization
encodeType(FreeAssAlgebra) = "FreeAssAlgebra"
@registerSerializationType(AbstractAlgebra.Generic.FreeAssAlgebra,true)



Oscar.serialize_with_id(FreeAssAlgebra) = true

function save_internal(s::Oscar.SerializerState, A::FreeAssAlgebra)
    d = Dict(
    :base_ring => Oscar.save_as_ref(s, base_ring(A)),
    :symbols => save_internal(s, symbols(A)),
    )
    return d
end


function load_internal(s::Oscar.DeserializerState, ::Type{<:FreeAssAlgebra}, dict::Dict)
    R = Oscar.load_unknown_type(s, dict[:base_ring])
    gens = Oscar.load_internal(s, Vector{Symbol}, dict[:symbols])
    return free_associative_algebra(R, gens)[1]
end


# Free associative algebra element serialization
encodeType(FreeAssAlgElem) = "FreeAssAlgElem"
@registerSerializationType(FreeAssAlgElem)





function save_internal(s::Oscar.SerializerState, f::FreeAssAlgElem)
    d = Dict(
    :parent => Oscar.save_as_ref(s, parent(f)),
    :coeffs => Oscar.save_internal(s, collect(coefficients(f))),
    :exps => save_internal(s, collect(exponent_words(f))),
    :length => save_internal(s, length(f)),
    )
    return d
end

function load_internal(s::Oscar.DeserializerState, ::Type{<:FreeAssAlgElem}, dict::Dict)
    parent = Oscar.load_ref(s, dict[:parent])
    coeff_type=elem_type(coefficient_ring(parent))
    coeffs = Oscar.load_internal(s, Vector{coeff_type}, dict[:coeffs]) 
    exps = Oscar.load_internal(s, Vector{Vector{Int}},  dict[:exps]) 
    length = Oscar.load_internal(s, Int, dict[:length])
    element = parent(coeffs, exps)
    return element
end


# Free associative algebra element serialization
encodeType(FreeAssAlgIdeal) = "FreeAssAlgIdeal"
@registerSerializationType(FreeAssAlgElem)





function save_internal(s::Oscar.SerializerState, f::FreeAssAlgElem)
    d = Dict(
    :parent => Oscar.save_as_ref(s, parent(f)),
    :coeffs => Oscar.save_internal(s, collect(coefficients(f))),
    :exps => save_internal(s, collect(exponent_words(f))),
    :length => save_internal(s, length(f)),
    )
    return d
end

function load_internal(s::Oscar.DeserializerState, ::Type{<:FreeAssAlgElem}, dict::Dict)
    parent = Oscar.load_ref(s, dict[:parent])
    coeff_type=elem_type(coefficient_ring(parent))
    coeffs = Oscar.load_internal(s, Vector{coeff_type}, dict[:coeffs]) 
    exps = Oscar.load_internal(s, Vector{Vector{Int}},  dict[:exps]) 
    length = Oscar.load_internal(s, Int, dict[:length])
    element = parent(coeffs, exps)
    return element
end









### Test stuff 
#
gens = ["x","y","z"]
A, g = FreeAssociativeAlgebra(Oscar.QQ, gens)
#
#elem_type(coefficient_ring(A))
#
#
#f = g[1]*g[2] + g[3]*g[1]
#coefficients(f)
e = collect(exponent_words(f))
#exponent_words(f) |> collect |> typeof
save("test",f)
#
#g = load("test")
#
#
#f.coeffs
#
#x = save_internal(Oscar.SerializerState(), base_ring(A))
#
#save("test",A)
#load("test")
#
#A1 = load("test")








