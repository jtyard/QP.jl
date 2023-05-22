# Some orders and their properties

using QP
using Oscar

include("su2k.jl")

# julia> is_principal(su2levelk(4-2).D)
# (true, 4)

is_principal(su2levelk(5-2).D)
# (true, 5)

# julia> is_principal(su2levelk(6-2).D)
# (true, 4)

# julia> is_principal(su2levelk(7-2).D)
# (true, 1)

# julia> is_principal(su2levelk(8-2).D)
# (true, 2)


