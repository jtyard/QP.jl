using Oscar, QP


F,z = cyclotomic_field(3)
mat(v) = matrix_space(F,3,3)([v[1] + v[2]*z, v[3] + v[4]*z, v[5] + v[6]*z,  v[7] + v[8]*z,v[9] + v[10]*z,v[11] + v[12]*z,v[13] + v[14]*z,v[15] + v[16]*z,v[17] + v[18]*z])

box(N) = Base.Iterators.product(repeat([-N:N],18)...)

I = matrix_space(F,3,3)([1,0,0,0,1,0,0,0,1])

U3 = [mat(v) for v in box(1) if mat(v)*dagger(mat(v))==I]
SU3 = [g for g in U3 if det(g) == 1]

length(U3), length(SU3)
#(1296, 216)



#length([mat(v) for v in box(2) if mat(v)*transpose(mat(v))==4*I])
# 8 so 