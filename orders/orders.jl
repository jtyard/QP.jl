# Some orders and their properties

using QP
using Oscar

include("su2k.jl")

k = 3
A = su2levelk(k)
OA = maximal_order(A)
OK = base_ring(OA)

D = discriminant(OA)

# julia> is_principal(su2levelkD(4-2))
# (true, 4)

# julia> is_principal(su2levelkD(5-2))
# (true, 5)

# julia> is_principal(su2levelkD(6-2))
# (true, 4)

# julia> is_principal(su2levelkD(7-2))
# (true, 1)

# julia> is_principal(su2levelkD(8-2))
# (true, 2)


