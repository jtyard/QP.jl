using Oscar, QP

N = 4

S = SicData(4)
F = S.F

P = projective_matrix_space(QQ,N)

I = Im(P) + Ihplus(P) + Ireal(P)



# Actually I think there are no real solutions but maybe a different symmetry is needed

#T = subscheme(P.P,I)
#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))

#rational_solutions(I)


#@time pd = primary_decomposition(I)