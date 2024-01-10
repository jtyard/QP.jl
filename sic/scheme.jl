using Oscar, QP

N = 3

P = projective_matrix_space(QQ,N)

X = gen(P)

laplacian(P,hplus(P,1,1))


Treal = subscheme(P.P,Im(P) + Ihplus(P) + Ireal(P))
T1 = subscheme(P.P,Im(P) + Ihplus(P))

Tall = subscheme(P.P,Ihminus(P) + Ihplus(P))

# dim(Tall)=1,5,7,13 for N=2,3,4,5 but N=6 takes a long time. 


function alldim(N)
    P = projective_matrix_space(QQ,N)
    dim(subscheme(P.P,Ihminus(P) + Ihplus(P)))
end

function alldimreal(N)
    P = projective_matrix_space(QQ,N)
    dim(subscheme(P.P,Ihminus(P) + Ihplus(P) + Ireal(P)))
end