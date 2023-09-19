###############
# Generalized Paulis and Heisenberg group
###############

export gpX, gpZ, heis, heis2, heiscocycle, heispairing, gauss_sum, weil_w0, weil_N, weil_T, weil_U

export heisAAZ, heisQ

export ABN, weil_overlaps, weil_ad, weil_zw

# generalized Pauli X 
function gpX(N::Int)
    X = zero_matrix(cyclotomic_field(N)[1],N,N)
    for i in 0:N-1
        X[1 + ((i+1) % N),1+i] = 1
    end
    X
end

# generalized Pauli Z
gpZ(N::Int) = diagonal_matrix([zetaN(N)^i for i in 0:N-1]...)

# Section of the single qudit Heisenberg group.  
function heis(j1::Int,j2::Int,N::Int)
    if iseven(N)
        C,z = cyclotomic_field(2N)
        (-z)^Int(-j1*j2)*map(C,gpZ(N)^(-j1)*gpX(N)^(-j2))
    else  
        twoinv = ZN(N)(2)^-1
        zetaN(N)^Int(-j1*j2*twoinv)*gpZ(N)^(-j1)*gpX(N)^(-j2)
    end
end

# It is related to the following section from
# eqn (8)  https://arxiv.org/abs/1209.1813 and eqn (6) of https://arxiv.org/abs/1604.06098
function heisAAZ(i::Int,j::Int,N::Int)
    C,z = isodd(N) ? cyclotomic_field(N) : cyclotomic_field(2*N);
    w = isodd(N) ? z : z 
    w^Int(i*j)*map(C,gpX(N)^Int(i)*gpZ(N)^Int(j))
end

# Section of multi-qudit heisenberg group
function heis(j::nmod_mat)
    N = Int(characteristic(base_ring(j)))
    n2 = length(j)
    if isodd(n2)  
        error("argument must have even length")
    end
    n = Int(n2//2)
    tensor([heis(Int(j[2i-1]),Int(j[2i]),N) for i in 1:n]...)
end

function heiscocycle(j::nmod_mat,k::nmod_mat) 
    N = Int(characteristic(base_ring(j))) 
    twoinv = ZN(N)(2)^-1
    w = zetaN(N)^Int(twoinv)  
    w^Int(j[1]*k[2] - j[2]*k[1])
end

function heispairing(j::nmod_mat,k::nmod_mat) 
    N = Int(characteristic(base_ring(j))) 
    zetaN(N)^Int((j[1]*k[2] - j[2]*k[1]))
end

# Alternatively we can always work modulo 2N
function heis2(j1::Int,j2::Int,N::Int)
    C,w = cyclotomic_field(2N)
    (-w)^(-j1*j2)*map(C,gpZ(N)^(-j1)*gpX(N)^(-j2))
end

function heis2(j::nmod_mat)
    N2 = Int(characteristic(base_ring(j)))
    n2 = length(j)
    
    @assert iseven(N2) "Characteristic must be even"
    @assert iseven(n2) "Argument must have even length"
    
    N = Int(N2//2)
    n = Int(n2//2)
    tensor([heis2(Int(j[2i-1]),Int(j[2i]),N) for i in 1:n]...)
end

# Quaternion section (N=2 only) for comparison
function heisQ(j::nmod_mat)
    N = Int(characteristic(base_ring(j)))
    @assert N==2 "Characteristic must equal 2"
    C,i = cyclotomic_field(4)
    #return ( i^(Int(j[1]*j[2])) ) * ( (C[0 i; i 0])^Int(j[1]) ) * ( (C[i 0; 0 -i])^Int(j[2]) )
    return  (C[0 i; i 0])^Int(j[1]) * (C[i 0; 0 -i])^Int(j[2])
end