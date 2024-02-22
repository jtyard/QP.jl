using Oscar, QP

S = SicData(5)

F = absolute_simple_field(S.F)[1]

C = cyclotomic_field(60)[1]

Ls = [L for L in subfields(F) if degree(L[1]) == 16]

L, L_to_F = Ls[findfirst(L -> is_isomorphic(L[1],C), Ls)]

L_to_C = is_isomorphic_with_map(L,C)[2]

FL = relative_simple_extension(F,L)[1]

OFL = maximal_order(FL)

DL = discriminant(OFL)

DC = L_to_C(DL)

p1, p2 = [p[1] for p in collect(factor(DC))]

print(DC == p1*p2)