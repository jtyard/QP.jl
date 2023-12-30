using Oscar, QP

#N = 2

S = projective_matrix_space(QQ,N)

X = gen(S)

PP = subscheme(S.P,Im(S)+Ihplus(S))

dim(PP)

