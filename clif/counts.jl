using Oscar, QP


F,i = cyclotomic_field(4)

mat(v) = matrix_space(F,2,2)([v[1] + v[2]*i, v[3] + v[4]*i, v[5] + v[6]*i,  v[7] + v[8]*i])

box(N) = Base.Iterators.product(repeat([-N:N],8)...)

I = matrix_space(F,2,2)([1,0,0,1])

U2 = [mat(v) for v in box(1) if mat(v)*dagger(mat(v))==I]

SU2 = [g for g in U2 if det(g) == 1]

O2 = [u for u in U2 if dagger(u) == transpose(u)]

SO2 = [u for u in O2 if det(u)==1]


length(U2), length(SU2), length(O2), length(SO2)
#(32, 8, 8, 4)
#(96, 24, 8, 4)


#length([mat(v) for v in box(2) if mat(v)*transpose(mat(v))==4*I])
# 8 so 