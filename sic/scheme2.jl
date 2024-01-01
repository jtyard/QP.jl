using Oscar, QP

N = 2

C12 = cyclotomic_field(12)[1]

P = projective_matrix_space(C12,N)

I = Im(P) + Ihplus(P)

#T = subscheme(P.P,I)
#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))

rational_solutions(I)


