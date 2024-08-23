

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
include("Poly/Polys.jl")

# A better approach to matrix polynomials?
include("Poly/MatPoly.jl")

# Cartesian products of polynomial rings 
include("Poly/ProdPoly.jl")

# Groups
include("Base/Groups.jl")

# Quadratic number fields and orders
include("Base/Quad.jl")

# Quaternion algebras and orders 
include("Quat/QuaternionAlgebras.jl")
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
