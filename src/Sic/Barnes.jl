# Barnes double gamma etc
using Oscar

setprecision(4096)
 
RRR() = ArbField(precision(BigFloat))
CCC() = AcbField(precision(BigFloat))
ii = onei(CCC())

RR = RRR()
CC(x) = CCC()(x)
CC(M::Array) = map(x->CC(x),M)

#import Nemo.erf
#import Nemo.exppii
#import Base.abs2
#import Base.*
#import Base./
# ∀ ⇒ ∪ ⋊ |+⟩|0⟩ + |-⟩|1⟩ ⟨⟩⊗ ψ ≈
#arc = Union{arb,acb}
#abs2(x::arc) = x*conj(x)

#erf(x::arb) = erf(CC(x))
#Erf(x) = erf(sqrt(RR(pi))*x)/2 # int_0^x e^{-pi u^2} du

#exppii(x::Float64) = exppii(CC(x))
#exppii(x::Int64) = exppii(CC(x))
#exppii(x::arb)     = exppii(CC(x))
#exp2pii(x)         = exppii(2*x)
#
#a::arc * M::Array = map(x->a*x, M)
#M::Array * a::arc = map(x->a*x, M)
#M::Array / a::arc = map(x->x/a, M)
#
#dot(u,v) = transpose(u)*v
#u⋅v = dot(u,v)

barnes