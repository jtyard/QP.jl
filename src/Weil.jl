###############
# Generalized Paulis and Heisenberg group
###############

export gpX, gpZ, heis, heiscocycle, heispairing, gauss_sum, weil_w0, weil_N, weil_T, weil

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


# Section of the Heisenberg group from AAZ & AFMY
# Do I need to change the signs/order?  What about for mod d and/or integer coordinates?
function heis(j::nmod_mat)
    N = Int(characteristic(base_ring(j)))
    if iseven(N)
        C,z = cyclotomic_field(2N)
        (-z)^Int(-j[1]*j[2])*map(C,gpZ(N)^Int(-j[1])*gpX(N)^Int(-j[2]))
    else  
        twoinv = ZN(N)(2)^-1
        zN(N)^Int(-j[1]*j[2]*twoinv)*gpZ(N)^Int(-j[1])*gpX(N)^Int(-j[2])
    end
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

###############
# Weil representation 
###############

function gauss_sum(a::nmod)
    N = Int(characteristic(parent(a)))
    z = cyclotomic_field(N)[2]
    a = Int(a)
    sum([z^(a*i^2) for i in 0:N-1])
end
    
# represent the Weyl element w0
function weil_w0(N::Int)
    F,z = cyclotomic_field(N)
    matrix(F, [[z^(-i*j) for j in 0:N-1] for i in 0:N-1])*(gauss_sum(ZN(N)(-2))*F(N)^(-1))
end

# represent the upper triangular unipotent matrices 
function weil_N(b::nmod)
    N = Int(characteristic(parent(b)))
    _,z = cyclotomic_field(N)
    b = Int(b)
    diagonal_matrix([z^(b * i^2 * invmod(2,N)) for i in 0:N-1]...)
end

# represent diagonal matrices
function weil_T(a::nmod)
    N = Int(characteristic(parent(a)))
    a = Int(a)
    M = zero_matrix(cyclotomic_field(N)[1],N,N)
    for i in 0:N-1
        M[1+i,1+a * i % N] = 1
    end
    jacobi_symbol(a,N)*M
end

# Accepts a matrix in SL(2,GF(p)), with p prime, and returns a p x p matrix over the pth cyclotomic field
function weil(g::nmod_mat)
    N = Int(characteristic(base_ring(g)))
    a = g[1,1]
    b = g[1,2]
    c = g[2,1]
    d = g[2,2]
    
    if c != 0 
        return weil_N(a*c^-1)*weil_w0(N)*weil_T(c)*weil_N(d*c^-1)
    else
        return weil_T(a)*weil_N(b*a^-1)
    end
end
