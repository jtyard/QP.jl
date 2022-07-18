# Groups

export eye

export PSL,PGL, abelian_2cocycle, abelian_2cocycle_mat
#export PGL, PSL, eye

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


