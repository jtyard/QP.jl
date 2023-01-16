module QP

using Oscar


# 
include("Base/Numbers.jl")
include("Base/Polys.jl")
include("Base/ElMats.jl")
include("Base/Groups.jl")
include("Base/Algebras.jl")

# Heisenberg groups for all Z/N and Weil for Z/p
include("Clif/Weil.jl")

# SIC-POVMs
include("Sic/Sic.jl")
# Quadratic number fields and orders
include("Sic/Quad.jl")
# Class fields

# Schur transform (eventually)
include("Schur/SU2.jl")


# Todo: Tensor categories
#include("Base/Cats.jl")
#include("Base/Types.jl")



end # module
