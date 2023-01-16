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

# Section of the Heisenberg group from AAZ & AFMY
# Do I need to change the signs/order?  What about for mod d and/or integer coordinates
function heis(j::nmod_mat)
    N = Int(characteristic(base_ring(j)))
    if iseven(N)
        C,z = cyclotomic_field(2N)
        (-z)^Int(-j[1]*j[2])*map(C,gpZ(N)^Int(-j[1])*gpX(N)^Int(-j[2]))
    else  
        twoinv = ZN(N)(2)^-1
        zetaN(N)^Int(-j[1]*j[2]*twoinv)*gpZ(N)^Int(-j[1])*gpX(N)^Int(-j[2])
    end
end

heis(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = heis(ZN(N)[Int(i) Int(j)])

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

# Let's work with 
function heis2(j::nmod_mat)
    N2 = Int(characteristic(base_ring(j)))
    @assert iseven(N2) "Characteristic must be even"
    N = Int(N2//2)
    C,w = cyclotomic_field(N2)
        (-w)^Int(j[1]*j[2])*map(C,gpX(N)^Int(j[1])*gpZ(N)^Int(j[2]))
end

# eqn (8) of AAZ
function heisAAZ(i,j,N)
    C,z = isodd(N) ? cyclotomic_field(N) : cyclotomic_field(2*N);
    w = isodd(N) ? z : z 
    w^Int(i*j)*map(C,gpX(N)^Int(i)*gpZ(N)^Int(j))
end


# Quaternion section (N=2 only)
function heisQ(j::nmod_mat)
    N = Int(characteristic(base_ring(j)))
    @assert N==2 "Characteristic must equal 2"
    C,i = cyclotomic_field(4)
    #return ( i^(Int(j[1]*j[2])) ) * ( (C[0 i; i 0])^Int(j[1]) ) * ( (C[i 0; 0 -i])^Int(j[2]) )
    return  (C[0 i; i 0])^Int(j[1]) * (C[i 0; 0 -i])^Int(j[2])
end

###############
# Unitary lifting of Weil representation of SL(2,Z/N) for prime N
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
function weil_U(g::nmod_mat)
    N = Int(characteristic(base_ring(g)))
    @assert is_prime(N)
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

weil_U(g::MatrixGroupElem{nmod, nmod_mat}) = weil_U(g.elm)


## 

function ABN(N; type = 0)
    X = gpX(N)
    Z = gpZ(N)
    O = 0*X
    if type == 0
        A = X
        B = Z
    elseif type == 1
        A = [O X; -X O]
        B = [O Z; -Z O]
    elseif type == 2 
        A = [X O; O X]
        B = [Z O; O Z]
    elseif type == 3
        A = [O X; X O]
        B = [O Z; Z O]
    end
    matrix_group([A,B])
end

###
# U(1)_k has anyons <zeta_N,-1>.
# For even k is a TQFT with anyon group Z/k
# For odd k is a spin TQFT with anyons Z/2k  
# The fermion k = ψ generates a Z/2 1-form symmetry 
###


# Given mxn U, return the mn x mn matrix for X ↦ UXU^-1 in the X11 ... X1n ... X2n ... Xmn basis.  

function adj(U)

end
    

###
# Adjoint action X -> U X U^-1 in Aut_0(B) specified by 
# (Z/N)^2 ⋊ SL(2,Z/N) (odd N) and (Z/N)^2 ⋊ SL(2,Z/2N) (even N)
###
function weil_ad_vec(g::nmod_mat) 
    n = Int(characteristic(base_ring(g)))
    if iseven(n)
        N = Int(n//2)
        error("Not implented")
    else
        N = n
        F,z_N = cyclotomic_field(N)
        U = zero_matrix(QQ,N^2,N^2)
        for j1 in 0:N-1
            for j2 in 0:N-1
                j = ZN(N)[j1; j2]
                k = g + j
                #println(j1,j2) 
                k1, k2 = Int(k[1]), Int(k[2])
                U[1 + k1*N + k2, 1 + j1*N + j2] = 1
            end
        end
    end
    U
end

function weil_ad_mat(g::nmod_mat) 
    n = Int(characteristic(base_ring(g)))
    if iseven(n)
        N = Int(n//2)
        error("Not implented")
    else
        N = n
        F,z_N = cyclotomic_field(N)
        U = zero_matrix(QQ,N^2,N^2)
        for j1 in 0:N-1
            for j2 in 0:N-1
                j = ZN(N)[j1; j2]
                k = g*j
                #println(j1,j2) 
                k1, k2 = Int(k[1]), Int(k[2])
                U[1 + k1*N + k2, 1 + j1*N + j2] = 1
            end
        end
    end
    U
end

function weil_ad(g::nmod_mat) 
    rc = nrows(g), ncols(g)
    if rc == (2,2) 
        return weil_ad_mat(g)
    elseif rc == (2,1) 
        return weil_ad_vec(g)
    elseif rc == (1,2)
        return weil_ad_vec(transpose(g))
    end
    error("Expects a matrix or a vector.")

end

weil_ad(g::MatrixGroupElem{nmod, nmod_mat}) = weil_ad(g.elm)

###
# Acts on the overlaps alpha_j(X) by (signed when N is even) permutation matrices.
# Specified by (Z/N)^2 ⋊ SL(2,Z/N) (odd N) and (Z/N)^2 ⋊ SL(2,Z/2N) (even N)
###

function weil_overlaps_vec(g::nmod_mat) 
    n = Int(characteristic(base_ring(g)))
    if iseven(n)
        N = Int(n//2)
        error("Not implented")
    else
        N = n
        F,z_N = cyclotomic_field(N)
        U = zero_matrix(QQ,N^2,N^2)
        for j1 in 0:N-1
            for j2 in 0:N-1
                j = ZN(N)[j1; j2]
                k = g + j
                #println(j1,j2) 
                k1, k2 = Int(k[1]), Int(k[2])
                U[1 + k1*N + k2, 1 + j1*N + j2] = 1
            end
        end
    end
    U
end



function weil_overlaps_mat(g::nmod_mat) 
    n = Int(characteristic(base_ring(g)))
    if iseven(n)
        N = Int(n//2)
        error("Not implented")
    else
        N = n
        F,z_N = cyclotomic_field(N)
        U = zero_matrix(QQ,N^2,N^2)
        for j1 in 0:N-1
            for j2 in 0:N-1
                j = ZN(N)[j1; j2]
                k = g*j
                #println(j1,j2) 
                k1, k2 = Int(k[1]), Int(k[2])
                U[1 + k1*N + k2, 1 + j1*N + j2] = 1
            end
        end
    end
    U
end


function weil_overlaps(g::nmod_mat) 
    rc = nrows(g), ncols(g)
    if rc == (2,2) 
        return weil_overlaps_mat(g)
    elseif rc == (2,1) 
        return weil_overlaps_vec(g)
    elseif rc == (1,2)
        return weil_overlaps_vec(transpose(g))
    end
    error("Expects a matrix or a vector.")

end

# Works with elements of SL(2,ZN(N)) too
weil_overlaps(g::MatrixGroupElem{nmod, nmod_mat}) = weil_overlaps(g.elm)

# G = SL(2,Z3); g1 = rand(G); g2 = rand(G); weil_ad(g1)*weil_ad(g2) == weil_ad(g1*g2)

# g = rand(G); j = rand(Z3^2); weil_ad(g)*weil_ad(j)*weil_ad(g^-1) == weil_ad(g*j)

# weil_U(g)*heis(j)*weil_U(g^-1) == heis(g*j)

