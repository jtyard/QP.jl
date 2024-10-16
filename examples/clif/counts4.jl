using Oscar, QP

n = 2
F,i = cyclotomic_field(4)

#mat(v) = matrix_space(F,2^n,2^n)([v[2*j-1] + i*v[2*j] for j in 1:4^n])


#box(N) = Base.Iterators.product(repeat([-N:N],2*4^n)...)

#I = identity_matrix(F,2^4)
#U2 = [mat(v) for v in box(1) if mat(v)*dagger(mat(v))==I]
#SU2 = [g for g in U2 if det(g) == 1]
#O2 = [u for u in U2 if dagger(u) == transpose(u)]
#SO2 = [u for u in O2 if det(u)==1]

#length(U2), length(SU2), length(O2), length(SO2)


O4 = automorphism_group(integer_lattice(;gram = identity_matrix(QQ,4)))

order(O4)
# 384

#order(automorphism_group(hermitian_lattice(F; gram = identity_matrix(F,2))))
#32

#order(automorphism_group(hermitian_lattice(F; gram = identity_matrix(F,4))))
#6144

#factor(6144)
#1 * 2^11 * 3