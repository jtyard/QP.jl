using Oscar, QP

N = 7

S = SicData(7)

P = S.P
X = gen(P)

I = Iminors(P,2) + Ihplus(P)  + Icyclic(P,Z7[2 0; 0 4])
J = Iminors(P,2) + Ihplus(P)  + Icyclic(P,Z7[3 0; 0 5])
Ir = Iminors(P,2) + Ihplus(P) + Ireal(P) + Icyclic(P,Z7[2 0; 0 4])
# Actually I think there are no real solutions but maybe a different symmetry is needed

T = subscheme(P,I)
#Tall = subscheme(P.P,Ihminus(S) + Ihplus(P))

#rational_solutions(I)

#mp = minimal_primes(I)d