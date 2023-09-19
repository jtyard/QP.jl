

module QP

using Oscar

# Scalars
include("Base/Numbers.jl")

# Elementary matrices and tensor products
include("Base/ElMats.jl")

# Polynomial rings
include("Base/Polys.jl")

# Groups
include("Base/Groups.jl")

# Quadratic number fields and orders
include("Base/Quad.jl")

# Quaternion algebras and orders 
include("Base/Quat.jl")

# Heisenberg groups for all Z/N and Weil for Z/p
include("Clif/Heis.jl")
include("Clif/Weil.jl")

# SIC-POVMs
include("Sic/Sic.jl")

# Abstract spin kets
include("Schur/SU2.jl")

end # module
