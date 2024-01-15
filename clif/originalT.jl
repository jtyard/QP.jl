using Oscar, QP


C12,z12 = cyclotomic_field(12)
s3 = sqrt(C12(-3))
z3 = z12^4
i = z12^3
mat(v) = matrix_space(C12,2,2)([v[1] + v[2]*z3, v[3] + v[4]*z3, v[5] + v[6]*z3,  v[7] + v[8]*z3])

box(N) = Base.Iterators.product(repeat([-N:N],8)...)

I = matrix_space(C12,2,2)([1,0,0,1])

#U2 = [mat(v) for v in box(2) if mat(v)*dagger(mat(v))==9*I]
#SU2 = [g for g in U2 if det(g) == 1]

#length(U2), length(SU2)
# (72,12) for N=1
# (,) for N=3

t12 = C12[3 1-i; 1+i -3]

x = C12[0 1; 1 0]
y = C12[0 -i;i 0]
z = C12[1 0; 0 -1]

u = (1 + i*(x+y+z))*C12(2)^-1