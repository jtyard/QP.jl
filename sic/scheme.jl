using Oscar, QP

N = 2

S = projective_matrix_space(QQ,N)

X = gen(S)

laplacian(S,hplus(S,1,1))


