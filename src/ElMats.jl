#################
# Elementary matrices
#################
export Eij, SWAP, Pplus, Pminus


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