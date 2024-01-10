using Oscar, QP

N = 2

C12 = cyclotomic_field(12)[1]

P12 = projective_matrix_space(C12,N)

I12 = Im(P12) + Ihplus(P12)

#T = subscheme(P.P,I)
#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))

#rational_solutions(I12)

P = projective_matrix_space(QQ,N)

I = Im(P) + Ihplus(P)

Ir = Ireal(P) + I

Iall = Ihplus(P) + Ihminus(P)

Iallr = Ihplus(P) + Ihminus(P) + Ireal(P)


#T = subscheme(P.P,I)
#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))
