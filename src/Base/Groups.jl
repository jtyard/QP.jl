# Groups
# Some useful / experimental routines for working with groups. 

using Oscar


export aut, inn

export eye, PSL,PGL, abelian_2cocycle, abelian_2cocycle_mat, affine_group 


import Oscar.aut
aut(G::Oscar.GAPGroup) = automorphism_group(G)
inn(G::Oscar.GAPGroup) = inner_automorphism_group(aut(G))[1]


#Given projective representation f : G -> matrices, return the corresponding normalized 2-cocycle
function abelian_2cocycle(G,f) 
    I = f(G(0))
    return (j,k) -> tr(f(j)*f(k)*f(j+k)^-1)//tr(f(G(0)))
end

function abelian_2cocycle_mat(G,f) 
    return matrix([abelian_2cocycle(G,f)(a,b) for a in G, b in G])
end

eye(G::MatrixGroup) = gens(G)[1]^0

function PGL(n,N) 
    A,u = unit_group(ZN(N))
    G = GL(n,ZN(N))
    H = sub([G([i==j ? u(a) : 0 for i = 1:n, j = 1:n]) for a in A]...)[1]
    quo(G,H)[1]
end

function PSL(n,N) 
    A,u = unit_group(ZN(N))
    B = [a for a in A if u(a)^n==1]
    G = SL(n,ZN(N))
    H = sub([G([i==j ? u(b) : 0 for i = 1:n, j = 1:n]) for b in B]...)[1]
    quo(G,H)[1]
end

function affine_group(G::MatrixGroup)
    n = degree(G)
    F = base_ring(G)
    linear_gens = [ [zero_matrix(F,1,n+1); zero_matrix(F,n,1) g.elm] for g in gens(G)]
    for i in 1:length(linear_gens)
        linear_gens[i][1,1] = F(1)
    end
    affine_gens = [identity_matrix(F,n+1) for i in 1:n]
    for i in 1:n
        affine_gens[i][i+1,1] = F(1)
    end
    return matrix_group(vcat(linear_gens,affine_gens))
end

