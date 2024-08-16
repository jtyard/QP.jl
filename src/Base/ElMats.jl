#################
# Elementary matrices
#################

#import Oscar.Hecke.⊗ as hecketensor

export Eij, SWAP, Pplus, Pminus, werner, bra, ket, tensor, tensorperm #, ⊗
export tr1, tr2, partial_trace, reduced_operator

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




    
function tr1(A::MatElem, d::Vector{Int})
    @assert length(d)==2
    m,n=d
    r, c = size(A)
    @assert r==c==m*n
    B = zero_matrix(base_ring(A),n,n)
    for i = 0:m-1
        B = B + A[n*i+1:n*i+n,n*i+1:n*i+n]
    end
    B
end


function tr2(A::MatElem, d::Vector{Int})
    @assert length(d)==2
    m,n = d
    r, c = size(A)
    @assert r==c==m*n
    B = zero_matrix(base_ring(A),m,m)
    for i = 0:m-1
        for j = 0:m-1
            B[i+1,j+1] = tr(A[n*i+1:n*i+n,n*j+1:n*j+n])
        end
    end
    B
end

#Helper functions for general partial trace
function indices_to_index(I::Vector{Int}, d::Vector{Int})
    n = length(d)
    @assert length(I) == n
    @assert all(0 ≤ I[a] ≤ d[a]-1 for a in 1:n)
    j = I[1]
    for a in 2:n
        j = d[n+2-a]*j + I[a] 
    end
    j + 1
end

function index_to_indices(i::Int, d::Vector{Int})
    n = length(d)
    @assert 1 ≤ i ≤ prod(d)
    J = Vector{Int}(undef, n)
    i -= 1
    for a in n:-1:1
        J[a] = i % d[a]
        i = div(i, d[a])
    end
    J
end

# Mother of all partial trace functions 
# Computes Tr_S i.e. the partial trace over the systems in S.  d contains the dimensions of the subsystems.
function partial_trace(A::MatElem, d::Vector{Int}, S::Vector{Int})
    n = length(d)
    @assert all(1 ≤ s ≤ n for s in S)
    @assert length(unique(S)) == length(S)
    
    total_dim = prod(d)
    @assert size(A) == (total_dim, total_dim)
    
    # Dimensions of the subsystem we're tracing out
    traced_dims = [d[s] for s in S]
    traced_dim = prod(traced_dims)

    # Dimensions of the subsystem we're keeping
    Sc = setdiff(1:n, S)
    kept_dims = [d[i] for i in Sc]
    kept_dim = prod(kept_dims)
    
    B = zero_matrix(base_ring(A), kept_dim, kept_dim)
    IK = [0 for s in 1:n]
    JK = [0 for s in 1:n]
    for i in 1:kept_dim, j in 1:kept_dim, k in 1:traced_dim
        I = index_to_indices(i,kept_dims)
        J = index_to_indices(j,kept_dims)
        K = index_to_indices(k,traced_dims)
        for a in 1:length(Sc)
            IK[Sc[a]] = I[a] 
            JK[Sc[a]] = J[a]
        end
        for a in 1:length(S)
            IK[S[a]] = K[a]
            JK[S[a]] = K[a]
        end
        B[i,j] += A[indices_to_index(IK,d), indices_to_index(JK,d)]
    end
    B
end

# d is a dimension vector, A is of prod d[i] x prod d[i] and s is a number in 1:length(d)
reduced_operator(A::MatElem, d::Vector{Int}, S::Vector{Int}) = partial_trace(A,d,setdiff(1:length(d),S))

# d is a dimension vector, A is of prod d[i] x prod d[i] and s is a number in 1:length(d)
reduced_operator(A::MatElem, d::Vector{Int}, s::Int) = reduced_operator(A,d,[s])

