#################
# Elementary matrices
#################

#import Oscar.Hecke.⊗ as hecketensor

export Eij, SWAP, Pplus, Pminus, werner, bra, ket, tensor, tensorperm #, ⊗
export partial_trace

# Matrix units
Eij(i::Int,j::Int,N::Int) = matrix(ZZ,[[a == (i % N) && b == (j % N) for b in 0:N-1] for a in 0:N-1])

# Access with an zzModRingElem vector like Eij(ZZmod(N)[i j])
Eij(ij::zzModMatrix) = Eij(Int(ij[1]), Int(ij[2]), Int(characteristic(base_ring(ij))))

# 2d array of matrix units 
Eij(N::Int) = [[Eij(i,j,N) for j in 0:N-1] for i in 0:N-1]

#################
# Unit vectors 
# Consider keeping these abstract like ket(Spin(1),Spin(0)), which 
# displays as |1,0⟩, and defining e.g. Vector{ZZ}(ket(...)) as below.
#############
ket(i::Int,N::Int) = matrix(ZZ,[[a==mod(i,N)] for a in 0:N-1])

ket(i::zzModRingElem) = ket(Int(i),Int(characteristic(parent(i))))

ket(a::zzModMatrix) = tensor([ket(a[i]) for i in 1:length(a)]...)

ket(a::FinGenAbGroupElem) = tensor([ket(Int(a[i]),Int(order(gens(parent(a))[i]))) for i in 1:length(a.coeff)]...)

bra(a...) = transpose(ket(a...))

# Tensor products

#tensor = kronecker_product

function Oscar.tensor(A...)
    n = length(A)
    out = A[n]
    for a in A[n-1:-1:1]
        out = kronecker_product(a,out)
    end
    out
end

#import Hecke.⊗
#Hecke.⊗=tensor

#A::MatElem ⊗ B::MatElem = tensor(A,B)
# Conflicts with Hecke.⊗ maybe import it directly?

# d is a vector of qudit dimensions, p is permutation of the same length
function tensorperm(d,p)
    e = [d[p(i)] for i in d]
    A = abelian_group(d)
    B = abelian_group(e)
    U = 0
    for a in A
        b = B([a[p(i)] for i in 1:length(a.coeff)])
        U = U + ket(b)*bra(a)
    end
    U
end

# Useful matrices on the tensor product

SWAP(N) = sum([kronecker_product(Eij(i,j,N),Eij(j,i,N)) for i in 1:N for j in 1:N])
Pplus(N) = (1//2)*(identity_matrix(QQ,N^2) + SWAP(N))
Pminus(N) = (1//2)*(identity_matrix(QQ,N^2) - SWAP(N))
werner(N,a) = (identity_matrix(QQ,N^2) + a*map(QQ,SWAP(N)))/(N^2 + a*N)


function partial_trace(A::MatElem, m,n)
    r, c = size(A)
    @assert r==c==m*n
    B = zero_matrix(base_ring(A),m,m)
    for i = 0:m-1
        for j = 0:m-1
            B[i+1,j+1] = tr(A[m*i+1:m*i+n,m*j+1:m*j+n])
        end
    end
    B
end

    
    