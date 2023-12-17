

module QP

using Oscar

nmod = zzModRingElem
nmod_mat = zzModMatrix
fmpz = ZZRingElem
fmpq = QQFieldElem
PolynomialRing = polynomial_ring

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
include("Sic/Fiducials.jl")

# Abstract spin kets
include("Schur/SU2.jl")

end # module
