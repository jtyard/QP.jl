using Oscar, QP

N = 5

S = SicData(5)

P = projective_matrix_space(QQ,5)

I = Im(P) + Ihplus(P)

T = subscheme(P.P,I)
#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))

@time dim(T)


