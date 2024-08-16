using Oscar, QP

N=3
R = MatPolyRing(QQ,N^3)
X = gen(R);

rho1 = partial_trace(X,[N,N,N],[2,3])
rho2 = partial_trace(X,[N,N,N],[1,3])
rho3 = partial_trace(X,[N,N,N],[1,2])

X1 = N*rho1 - rho1^0*trace(rho1)
X2 = N*rho2 - rho2^0*trace(rho2)
X3 = N*rho3 - rho3^0*trace(rho3)

@time I1 = ideal(R.Sgr,vec(collect(X1)));
@time I2 = ideal(R.Sgr,vec(collect(X2)));
@time I3 = ideal(R.Sgr,vec(collect(X3)));
@time Im = Iminors(R,2)
I = I1 + I2 + I3 + Im

@time dim(proj(I))