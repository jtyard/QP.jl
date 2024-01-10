using Oscar, QP

N = 7

S = SicData(7)

P = S.P

I = Im(P) + Ihplus(P)  + It(P,Z7[2 0; 0 4])
J = Im(P) + Ihplus(P)  + It(P,Z7[3 0; 0 5])
Ir = Im(P) + Ihplus(P) + Ireal(P) + It(P,Z7[2 0; 0 4])
# Actually I think there are no real solutions but maybe a different symmetry is needed

T = subscheme(P.P,I)
#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))

#rational_solutions(I)

#mp = minimal_primes(I)d