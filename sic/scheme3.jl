using Oscar, QP

N = 3

P = projective_matrix_space(rationals_as_number_field()[1],N)
R = P.RP

f = dwork_polynomial(P,0,0)

I = Im(P) + Ihplus(P) #+ Ireal(P)

II = I + ideal(R,f)

T = subscheme(P.P,I)
TT = subscheme(T,f)



#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))

#rational_solutions(I)


#@time pd = primary_decomposition(I)

dim(T), dim(subscheme(T,f)), is_saturated(I), is_saturated(II), is_radical(I), is_radical(II)
#(1, 0, true, true, false, false)