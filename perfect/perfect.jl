using Oscar

function tr2(A::MatElem, m,n)
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

function tr1(A::MatElem, m,n)
    r, c = size(A)
    @assert r==c==m*n
    B = zero_matrix(base_ring(A),n,n)
    for i = 0:m-1
        B = B + A[m*i+1:m*i+n,m*i+1:m*i+n]
    end
    B
end

P = projective_space(QQ,["z00", "z01", "z10", "z11", "w00", "w01", "w10", "w11"])

R = homogeneous_coordinate_ring(P)

(z00, z01, z10, z11, w00, w01, w10, w11) = gens(R)

psi = matrix([z00;z01;z10;z11])
psibar = matrix([w00;w01;w10;w11])

Psi = psi*transpose(psibar)

r = dot(psi,psibar)

rho1 = tr2(Psi,2,2)
rho2 = tr1(Psi,2,2)
eqns1 = collect(2*rho1 - r*identity_matrix(R,2))[:]
eqns2 = collect(2*rho2 - r*identity_matrix(R,2))[:]
#display(eqns)

S1 = subscheme(P,eqns1)
S2 = subscheme(P,eqns2)

I1 = defining_ideal(S1)
I2 = defining_ideal(S2)

S = subscheme(P,I1 + I2)

# is it is reduced
println(is_reduced(S))

# of dimension 4
println(dim(S))  

# and irreducible
println(is_irreducible(S))

#On the other hand each of S1 and S2 has 2 components each of dimension 4: 

# We can compute the dimensions of the components
J1 = [p[2] for p in primary_decomposition(I1,algorithm=:SY)]
J2 = [p[2] for p in primary_decomposition(I2,algorithm=:SY)]

println([dim(subscheme(P,p)) for p in J1])
println([dim(subscheme(P,p)) for p in J2])

