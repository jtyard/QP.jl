using Oscar, QP

m = 2
N = 4*m^2+3

S = SicData(N)

P = S.P
X = gen(P)
f = tr(X) - 1


r = 2*m^2 + 1
a = filter(x -> myorder(x) == r,collect(ZN(N)))[1]
t = ZN(N)[a 0; 0 a^-1]

I = Iminors(P,2) + Ihplus(P)  + Icyclic(P,t)
Ir = Iminors(P,2) + Ihplus(P) + Ireal(P) + Icyclic(P,t)
# Actually I think there are no real solutions but maybe a different symmetry is needed

T = subscheme(P,Ir)

A = affine_slice(T,f)

vdim(A)
