
#include("./utils.jl")
#include("./multiSetMatroid.jl")
#include("./save.jl")

export is_commutative,
    magic_unitary

@doc raw"""
magic_unitary(A::FreeAssociativeAlgebra)

This function takes a FreeAssociativeAlgebra and returns the magic unitary matrix associated with it.

# Examples

```julia
using Oscar
I = quantum_symmetric_group(3)
A = base_ring(I)
magic_unitary(A)
```
"""
function magic_unitary(A::FreeAssociativeAlgebra)
    n = Int(floor(sqrt(ngens(A))))
    @assert n^2 == ngens(A) "The algebra is not square"
    u = permutedims(reshape(gens(A),(n,n)),[2,1])
end

magic_unitary(I::Oscar.FreeAssociativeAlgebraIdeal) = magic_unitary(base_ring(I))


@doc raw"""
# Examples

```jldoctest
using Oscar
I = quantum_symmetric_group(4)
A = base_ring(I)

gb = groebner_basis(I)

save("test.test",collect(gb))
load("test.test")

magic_unitary(A)
QuantumAutomorphismGroups.is_commutative(I)

# output

false

```

"""
function is_commutative(I::Oscar.FreeAssociativeAlgebraIdeal; deg_bound::Int = -1, detailed::Bool = false)
  A = base_ring(I)
  u = magic_unitary(A)
  n = Int(floor(sqrt(ngens(A))))  
  
  commutator_matrix = fill(true,n,n)
  time_matrix = zeros(n,n)

  for i in 1:n, j in i+1:n
    commutator_matrix[i,j], time_matrix[i,j] = @timed ideal_membership(u[i,j]*u[j,i] - u[j,i]*u[i,j],I,deg_bound)
  end
 
  if detailed
    return commutator_matrix, time_matrix
  end
  
  return all(commutator_matrix)
end
