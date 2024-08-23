###############
# Multivariate polynomial rings of matrix entries
###############

using Oscar
using Unicode

import Oscar.gen, Oscar.matrix, Oscar.minors, Oscar.irrelevant_ideal, Oscar.saturation, Oscar.hom

export matrix_polynomial_ring, MatPolyRing, MatPolyIdeal


export Iminors, Ihplus, Ihminus, Itorus, dwork_modulus, dwork_polynomial

export MatPolyMap

include("MatPolyRing.jl")
include("MatPolyMap.jl")
#include("MatPolyIdeal.jl")


#######################################################################
#######################################################################
# Useful polynomials 
#######################################################################
#######################################################################


function fij(S::MatPolyRing,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    X(i,j) = gen(S,i,j)
    sum([X(a,a+j1)*X(a+j1+j2,a+j2) for a=0:S.N-1])
end

function fij(S::MatPolyRing,j::zzModMatrix) 
    fij(S,j[1],j[2])
end

function fplus(S::MatPolyRing,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    (1//2)*(fij(S,j1,j2) + fij(S,j2,j1))
end

function fplus(S::MatPolyRing,j::zzModMatrix) 
    fplus(S,j[1],j[2])
end

function fminus(S::MatPolyRing,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    (1//2)*(fij(S,j1,j2) - fij(S,j2,j1))
end

function fminus(S::MatPolyRing,j::zzModMatrix) 
    fminus(S,j[1],j[2])
end

## Projecting onto harmonic subspace

function hij(S::MatPolyRing,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    X = gen(S) 
    c = (1//2)*(tr(X^2) + tr(X)^2)*(1//(S.N+1))
    fij(S,j1,j2) - S.Sgr((Int(j1) == 0 ? 1 : 0) + (Int(j2) == 0 ? 0 : 1))*c
end

function hij(S::MatPolyRing,j::zzModMatrix) 
    hij(S,j[1],j[2])
end

function hplus(S::MatPolyRing,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem})
    X = gen(S) 
    c = (1//2)*(tr(X^2) + tr(X)^2)*(1//(S.N+1))
    fplus(S,j1,j2) - ((Int(j1) == 0 ? 1 : 0) + (Int(j2) == 0 ? 1 : 0))*c
end

function hplus(S::MatPolyRing,j::zzModMatrix) 
    hplus(S,j[1],j[2])
end

function hminus(S::MatPolyRing,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    X = gen(S) 
    c = (1//2)*(tr(X^2) - tr(X)^2)*(1//(S.N-1))
    fminus(S,j1,j2) - ((Int(j1) == 0) || (Int(j2) == 0) ? 1 : 0)*c
end

function hminus(S::MatPolyRing,j::zzModMatrix) 
    hminus(S,j[1],j[2])
end


## Laplacian 
# Should eventually replace Laplacian and laplacian from Polys.jl
function laplacian(S::MatPolyRing,f)
    sum([derivative(f,gen(S,a,a)) for a in 0:S.N-1])
end

## Ideals
## To make life easier for now they all live in the graded ring S.Sgr
## i.e. so they are not actually of type MatPolyIdeal
## Might need a MatDecPolyIdeal to do it properly...

function Ihplus(S::MatPolyRing)
    ideal(S.Sgr,[hplus(S,j1,j2) for j1 in 0:S.N-1 for j2 in 0:j1])
end

Ihminus(S::MatPolyRing) = ideal(S.Sgr,[hminus(S,i,j) for i in 0:S.N-1 for j in 0:i if i != j])

function Iminors(S::MatPolyRing,k) 
    ideal(S,minors(gen(S),k))
end

Ireal(S::MatPolyRing) = ideal(S.Sgr,[gen(S,i,j) - gen(S,j,i) for i in 0:S.N-1 for j in 0:i if i != j])

function Itorus(S::MatPolyRing,t::zzModMatrix) 
    N = S.N
    ideal(S.Sgr,[gen(S,j1,j2) - gen(S,Int((ZN(N)[j1 j2]*t)[1]),Int((ZN(N)[j1 j2]*t)[2])) for j1 in 0:S.N-1 for j2 in 0:S.N-1])
end

# Ideal generated by action of cyclic group <t> on generators.  
function Icyclic(S::MatPolyRing,t::zzModMatrix) 
    N = S.N
    ideal(S.Sgr,[gen(S,ZN(N)[j1 j2]) - gen(S,ZN(N)[j1 j2]*t) for j1 in 0:S.N-1 for j2 in 0:S.N-1])
end

#function Oscar.irrelevant_ideal(R::MPolyRing)
#    ideal(gens(R))
#end

Ic(R::MatPolyRing) = ideal(R.Sgr,[gen(R,j1,j2) - gen(R,j2,j1) for j1 in 0:R.N-1 for j2 in 0:R.N-1 if j1 < j2] )
Icc(R::MatPolyRing) = ideal(R.Sgr, [gen(R,j1,j2) - gen(R,-j2,-j1) for j1 in 0:R.N-1 for j2 in 0:R.N-1] )
Itr1(R::MatPolyRing) = ideal(R,[tr(gen(R)) - 1])
Itr0(R::MatPolyRing) = ideal(R.Sgr,[tr(gen(R))])
IT(R::MatPolyRing, a::Int) = ideal(R.Sgr, [gen(R)(j1,j2) - gen(R,a*j1,a*j2) for j1 in 0:R.N-1 for j2 in 0:R.N-1] )



function irrelevant_ideal(S::MatPolyRing)
    irrelevant_ideal(S.Sgr)
end


function dwork_polynomial(P::MatPolyRing,mu,mubar)
    X = gen(P)
    N = P.N
    
    ss = sum(X[i,j]^N for i in 1:N for j in 1:N) # (sum z_i^N)(sum w_i^N)
    pp = prod(X[i,i] for i in 1:N) # (prod z_i^N)(prod w_i^N)
    sp = sum(prod(X[i,j] for j in 1:N) for i in 1:N) # (sum z_i^N)(prod w_i^N)
    ps = sum(prod(X[j,i] for j in 1:N) for i in 1:N) # (prod z_i^N)(sum w_i^N)
    
    ss - N*(ps*mu + sp*mubar) + N^2*mu*mubar*pp
end