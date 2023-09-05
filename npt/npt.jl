# Can SICs find NPT bound entanglement?

using Oscar, QP

rho = werner(4,-1//2)

S4 = SicData(4)
F = S4.F 
c = S4.c  
S = [transpose(P[1,:]) for P in sic(4)]

U  = tensorperm([2,2,2,2],@perm (1, 3,2,4))

v = rand(S) 
w = rand(S)

p = tensor(v,v,v,v) + tensor(w,w,w,w)

a = (dagger(p,c)*tensor(rho,rho)*p)[1]

# next steps: pick a standard embedding of F and a root of unity in SicData
[e(a) for e in complex_embeddings(F)]

