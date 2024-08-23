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

#d = 2
P = projective_space(QQ,2*d^2-1)

R = homogeneous_coordinate_ring(P)

psi = matrix(gens(R)[1:d^2])
psibar = matrix(gens(R)[d^2+1:2*d^2])

Psi = psi*transpose(psibar)

r = dot(psi,psibar)

rho1 = tr2(Psi,d,d)
rho2 = tr1(Psi,d,d)
eqns1 = collect(d*rho1 - r*identity_matrix(R,d))[:]
eqns2 = collect(d*rho2 - r*identity_matrix(R,d))[:]
#display(eqns)

S = subscheme(P,vcat(eqns1,eqns2))
dim(S)
#S2 = subscheme(P,eqns2)
#
#I1 = defining_ideal(S1)
#I2 = defining_ideal(S2)

#S = subscheme(P,I1 + I2)

#is it is reduced
#println(is_reduced(S))

# of dimension 4
#println(dim(S))  

# and irreducible
#println(is_irreducible(S))

#On the other hand each of S1 and S2 has 2 components each of dimension 4: 

# We can compute the dimensions of the components
#J1 = [p[2] for p in primary_decomposition(I1,algorithm=:SY)]
#J2 = [p[2] for p in primary_decomposition(I2,algorithm=:SY)]

#println([dim(subscheme(P,p)) for p in J1])
#println([dim(subscheme(P,p)) for p in J2])

# The second components coincide. 
#J1[2] == J2[2]