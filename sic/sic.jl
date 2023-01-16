# Eventually put the relevant structures into SicData.  Seems to kind of work..

using Oscar
using QP


N = 4;

R = QQX(N,graded=true)
Zd = ZN(N)
Z2d = ZN(2*N)
I = Ih(N,graded=true) + Im(N,graded=true) 

S = ProjectiveScheme(R,I)
aff = affine_cone(S)
println(dim(aff)-1)

G = SL(2,ZN(N))
V = FreeModule(ZN(N),2)

g = rand(G)
v = rand(V)