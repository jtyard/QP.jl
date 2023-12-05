using Oscar, QP


F,i = cyclotomic_field(4)
mat(v) = matrix_space(F,2,2)([v[1] + v[2]*i v[3] + v[4]*i; v[5] + v[6]*i  v[7] + v[8]*i])

box(N) = Base.Iterators.product(repeat([-N:N],8)...)

[mat(v) for v in box(1) if mat(v)]