module QP

using Oscar
using Caching

export ZZN, zetaN, qint 

export Eij, SWAP, Pplus, Pminus




# abstract type QuantumSystem end




#################
# Ring   Element type
# ZZ     fmpz 
# QQ     fmpq
# ZN(N)  nmod
#
# Note that characteristic(base_ring(___)) works for GF(p)[ ] and ResidueClassRing(ZZ,4)
#################

Base.Int(a::nmod) = Int(ZZ(a))

ZZN(N) = ResidueRing(ZZ,N) 

zetaN(N) = cyclotomic_field(N)[2]

# However should rewrite this to live in the real cyclotomic field (also which version???)
qint(n,m) = sum([zetaN(2*n)^(m-1 - 2*i) for i in 0:m-1])


#################
# Matrix units and projectors
#################

# Matrix units
Eij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = matrix(ZZ,[[a == (Int(i) % N) && b == (Int(j) % N) for b in 0:N-1] for a in 0:N-1])

# Access with an nmod vector like Eij(ZZmod(N)[i j])
Eij(ij::nmod_mat) = Eij(ij[1], ij[2], Int(characteristic(base_ring(ij))))

# 2d array of matrix units 
Eij(N::Int) = [[Eij(i,j,N) for j in 0:N-1] for i in 0:N-1]

# Useful matrices on the tensor product

SWAP(N) = sum([kronecker_product(Eij(i,j,N),Eij(j,i,N)) for i in 1:N for j in 1:N])
Pplus(N) = (1//2)*(identity_matrix(QQ,N^2) + SWAP(N))
Pminus(N) = (1//2)*(identity_matrix(QQ,N^2) - SWAP(N))


include("Vars.jl")
include("Weil.jl")
include("Algebras.jl")


end # module
