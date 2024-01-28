using Oscar, QP

F,w = cyclotomic_field(3)

phi = F[0; 1; -1]

adder = sum([ket(Z3[a+b a-b])*bra(Z3[a b]) for a in Z3 for b in Z3])

I = identity_matrix(F,3)

