

module QP

using Oscar

#nmod = zzModRingElem
#nmod_mat = zzModMatrix
#fmpz = ZZRingElem
#fmpq = QQFieldElem
#PolynomialRing = polynomial_ring

# Scalars
include("Base/Numbers.jl")

# Elementary matrices and tensor products
include("Base/ElMats.jl")

# Polynomial rings (to be removed/replaced by ambient schemes)
include("Base/Polys.jl")

# Projective scheme of matrices over arbitrary base. Include harmonics and minors here.
include("Base/Schemes.jl")

# Groups
include("Base/Groups.jl")

# Quadratic number fields and orders
include("Base/Quad.jl")

# Quaternion algebras and orders 
include("Quat/Orders.jl")
include("Quat/ClassNumber.jl")

# Heisenberg groups for all Z/N and Weil for Z/p
include("Clif/Heis.jl")
include("Clif/Weil.jl")

# SIC-POVMs
include("Sic/Sic.jl")
include("Sic/Fiducials.jl")

# Abstract spin kets
include("Schur/SU2.jl")

end # module
