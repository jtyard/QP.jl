###############
# Generalized Paulis and Heisenberg group
###############

export gauss_sum, weil_w0, weil_N, weil_T, weil_U

export weil_overlaps, weil_ad, weil_zw

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
    @assert is_odd(N)
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


###
# U(1)_k Chern-Simons theory has anyons <zeta_N,-1>.
# For odd k is a 2+1D spin TQFT with anyons Z/2k and k = psi is the fermion.
# For even k is a 2+1D TQFT with anyon group Z/k.  There is also a spin-TQFT with transparent fermion and anyon group Z/k x Z/2. 
# In each case the fermion generates a Z/2 1-form symmetry that can be condensed to give a super TQFT (via gauging the fermion parity)
###


# Given mxn U, return the mn x mn matrix for X ↦ UXU^-1 in the X11 ... X1n ... X2n ... Xmn basis.  

function adj(U)
    error("Not implemented")
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

