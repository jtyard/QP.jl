module QP

using Oscar
using Caching

abstract type QuantumSystem end

# 
include("Base/Numbers.jl")
include("Base/Vars.jl")
include("Base/ElMats.jl")
include("Base/Groups.jl")
include("Base/Algebras.jl")

# Todo: Tensor categories
include("Base/Cats.jl")

# Heisenberg groups for all Z/N and Weil for Z/p
include("Clif/Weil.jl")

include("Sic/Sic.jl")

# Schur transform (someday)
include("Schur/Schur.jl")


end # module
