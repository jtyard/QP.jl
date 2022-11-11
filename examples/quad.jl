using Oscar

F,s = quadratic_field(-5)
OF = maximal_order(F)
s = OF(s)
C,ph = class_group(OF)

I3 = ideal(OF,OF(3))
P1 = collect(keys(factor(I3)))[1]
P2 = collect(keys(factor(I3)))[2]
I3 == P1*P2