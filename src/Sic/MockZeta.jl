# Will implement Kopp's "mock zeta functions" = completed indefinite zeta functions 
# Some of the first julia code I wrote, gave up on getting it to work but maybe look again? 
using Nemo

setprecision(4096)
 
RRR() = ArbField(precision(BigFloat))
CCC() = AcbField(precision(BigFloat))
ii = onei(CCC())

RR = RRR()
CC(x) = CCC()(x)
CC(M::Array) = map(x->CC(x),M)

import Nemo.erf
import Nemo.exppii
import Base.abs2
import Base.*
import Base./
# ∀ ⇒ ∪ ⋊ |+⟩|0⟩ + |-⟩|1⟩ ⟨⟩⊗ ψ ≈
arc = Union{arb,acb}
abs2(x::arc) = x*conj(x)

erf(x::arb) = erf(CC(x))
Erf(x) = erf(sqrt(RR(pi))*x)/2 # int_0^x e^{-pi u^2} du

exppii(x::Float64) = exppii(CC(x))
exppii(x::Int64) = exppii(CC(x))
exppii(x::arb)     = exppii(CC(x))
exp2pii(x)         = exppii(2*x)

a::arc * M::Array = map(x->a*x, M)
M::Array * a::arc = map(x->a*x, M)
M::Array / a::arc = map(x->x/a, M)

dot(u,v) = transpose(u)*v
u⋅v = dot(u,v)

# Indefinite theta function for g = 2 following Defn 4.10 of https://arxiv.org/abs/1912.12364
# Assumes Om is 2x2 symmetric and indefinite, real c = [c1,c2] with  c^T imag(Om) c < 0

Erfn(z,Om,c1,c2,n) = Erf(c2⋅(n+imag(z))/sqrt(-(c2⋅imag(Om)*c2)/2)) - Erf(c1⋅(n+imag(z))/sqrt(-(c1⋅imag(Om)*c1)/2)) 
Expn(z,Om,c1,c2,n) = exp2pii((n⋅Om*n)/2 + n⋅z)

function IndefiniteTheta(z,Om,c1,c2,N=10)
    z,Om,c1,c2 = CC(z),CC(Om),CC(c1),CC(c2)
    C1 = c1/sqrt(-(c1⋅imag(Om)*c1)/2)
    C2 = c2/sqrt(-(c2⋅imag(Om)*c2)/2)

    out = CC(0)
    for n1 in -N:N 
        for n2 in -N:N 
            n = [n1, n2] 
            out += ( Erf(C2⋅imag(Om*n + z)) - Erf(C1⋅imag(Om*n + z)) ) * exp2pii((n⋅Om*n)/2 + n⋅z)
        end
    end
    out
end

# Defn 4.15 
# p,q in R^2, Om in C^{2x2} symmetric with imag(Om) of signature (1,1), c1,c2 in R^2 of negative norm
IndefiniteThetaNull(p,q,Om,c1,c2,N=10)=exp2pii(q⋅Om*q/2+p⋅q)*IndefiniteTheta(p+Om*q,Om,c1,c2,N)

function XZ(j,d::Int) 
    out = fill(CC(0),d,d)
    for a in 0:d-1
        out[(a + j[1]) % d + 1,a + 1] = exp2pii(1/d)^(a*j[2])
    end
    out
end

function fj(z,j1,j2)
    d = length(z)
    zwwz = sum([z[1+a]*conj(z[1+(a+j1)%d])*conj(z[1+(a+j2)%d])*z[1+(a+j1+j2)%d] for a in 0:d-1])
    zwwz - ( (((j1 % d) == 0) ? 1 : 0) + (((j2 % d) == 0) ? 1 : 0) )*(conj(z)⋅z)^2/(d+1)
end

function checks(z)
    d = length(z)
    [[fj(z,j1,j2) for j2 in 0:d-1] for j1 in 0:d-1]/(conj(z)⋅z)^2
end

function check(z)
    d = length(z)
    sum([abs2(fj(z,j1,j2)) for j1 in 0:d-1 for j2 in 0:d-1])/real((conj(z)⋅z)^2)
end

# Test case from Kopp
Om = ii*[2  0 ; 0  (-6)] 
c1 = [0,1]
P = [2  3 ; 1  2]
c2 = P*c1
C1 = c1/sqrt(-(c1⋅imag(Om)*c1)/2)
C2 = c2/sqrt(-(c2⋅imag(Om)*c2)/2)
z = [0,0]
p0 = [0,0]
q0 = [1/5,0]
M = [2  0 ; 0  (-6)]

IndefiniteTheta(Om*q0,Om,c1,c2)

clam(lam) = (1-lam)*c1 + lam*c2
Alam(lam) = M + M*clam(lam)*transpose(clam(lam))*M/(-transpose(clam(lam))*M*clam(lam)/2)
Clam(lam) = clam(lam)/sqrt(-(clam(lam)⋅M*clam(lam))/2)