module QP

using Oscar

using Caching

export Oscar

export ZZmod

export Eij

export SWAP, Pplus, Pminus

export MatrixPolynomialRing, VariableMatrix

export QuaternionAlgebra

export qint

abstract type QuantumSystem end




#################
# Ring   Element type
# ZZ     fmpz 
# QQ     fmpq
# ZN(N)  nmod
#
# Note that characteristic(base_ring(___)) works for GF(p)[ ] and ResidueClassRing(ZZ,4)
#################

Base.Int(a::nmod) = Int(ZZ(a))

ZZmod(N) = ResidueRing(ZZ,N) 

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


###############
# Matrix polynomials
###############

function MatrixPolynomialRing(F,N::Int,X::Union{AbstractString, Char, Symbol} = "X")
    PolynomialRing(F,[string(X,"_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1],cached=true)[1]
end

function VariableMatrix(F,N::Int,X::Union{AbstractString, Char, Symbol} = "X")
    R = MatrixPolynomialRing(F,N,X)
    matrix(R,N,N,gens(R))
end


#########
# Quantum integers and things
#########

qint(n,m) = sum([CyclotomicField(2*n)[2]^(m-1 - 2*i) for i in 0:m-1])


#################
# Algebras and orders
#################

# Define a wrapper function that coerces the generators automatically 
QuaternionAlgebra(K,a,b) = Hecke.AlgQuat(K,K(a),K(b))



end
