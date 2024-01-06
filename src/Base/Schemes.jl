

using Oscar

import Oscar.gen, Oscar.matrix, Oscar.minors, Oscar.irrelevant_ideal, Oscar.saturation

export ProjectiveMatrixSpace, projective_matrix_space, gen, matrix

export fij, fplus, fminus, hij, hplus, hminus

export laplacian

export Iminors, Ihplus, Ihminus, Ireal, It, irrelevant_ideal, saturation, is_saturated


struct ProjectiveMatrixSpace
    P::ProjectiveScheme
    N::Int
    X::Union{AbstractString, Char, Symbol}
    RP::MPolyDecRing
    A::Spec 
    RA::MPolyRing
    ph::Oscar.MPolyAnyMap

    
    function ProjectiveMatrixSpace(F::AbstractAlgebra.Ring,N::Int;X::Union{AbstractString, Char, Symbol} = "X")
        P = projective_space(F,[string(X,"_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1])
        A,ph = affine_cone(P)
        new(P,N,X,homogeneous_coordinate_ring(P),A,coordinate_ring(A),ph)
    end
end

function projective_matrix_space(F::AbstractAlgebra.Ring,N::Int;X::Union{AbstractString, Char, Symbol} = "X")
    ProjectiveMatrixSpace(F,N,X=X)
end


matrix(S::ProjectiveMatrixSpace) = matrix(S.RP,S.N,S.N,gens(homogeneous_coordinate_ring(S.P)))

gen(S::ProjectiveMatrixSpace) = matrix(S)

function gen(S::ProjectiveMatrixSpace,i::Union{Int,zzModRingElem},j::Union{Int,zzModRingElem})
    #@assert divides(Int(S.N),Int(characteristic(base_ring(j))))[1]
    matrix(S)[Int(ZN(S.N)(i))+1,Int(ZN(S.N)(j))+1]
end

function gen(S::ProjectiveMatrixSpace,j::zzModMatrix)
    #@assert divides(Int(S.N),Int(characteristic(base_ring(j))))[1]
    gen(S,j[1],j[2])
end

function fij(S::ProjectiveMatrixSpace,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    X(i,j) = gen(S,i,j)
    sum([X(a,a+j1)*X(a+j1+j2,a+j2) for a=0:S.N-1])
end

function fij(S::ProjectiveMatrixSpace,j::zzModMatrix) 
    fij(S,j[1],j[2])
end

function fplus(S::ProjectiveMatrixSpace,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    (1//2)*(fij(S,j1,j2) + fij(S,j2,j1))
end

function fplus(S::ProjectiveMatrixSpace,j::zzModMatrix) 
    fplus(S,j[1],j[2])
end

function fminus(S::ProjectiveMatrixSpace,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    (1//2)*(fij(S,j1,j2) - fij(S,j2,j1))
end

function fminus(S::ProjectiveMatrixSpace,j::zzModMatrix) 
    fminus(S,j[1],j[2])
end

## Projecting onto harmonic subspace

function hij(S::ProjectiveMatrixSpace,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    X = gen(S) 
    c = (1//2)*(tr(X^2) + tr(X)^2)*(1//(S.N+1))
    fij(S,j1,j2) - S.RP((Int(j1) == 0 ? 1 : 0) + (Int(j2) == 0 ? 0 : 1))*c
end

function hij(S::ProjectiveMatrixSpace,j::zzModMatrix) 
    hij(S,j[1],j[2])
end

function hplus(S::ProjectiveMatrixSpace,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem})
    X = gen(S) 
    c = (1//2)*(tr(X^2) + tr(X)^2)*(1//(S.N+1))
    fplus(S,j1,j2) - ((Int(j1) == 0 ? 1 : 0) + (Int(j2) == 0 ? 1 : 0))*c
end

function hplus(S::ProjectiveMatrixSpace,j::zzModMatrix) 
    hplus(S,j[1],j[2])
end

function hminus(S::ProjectiveMatrixSpace,j1::Union{Int,zzModRingElem},j2::Union{Int,zzModRingElem}) 
    X = gen(S) 
    c = (1//2)*(tr(X^2) - tr(X)^2)*(1//(S.N-1))
    fminus(S,j1,j2) - ((Int(j1) == 0) || (Int(j2) == 0) ? 1 : 0)*c
end

function hminus(S::ProjectiveMatrixSpace,j::zzModMatrix) 
    hminus(S,j[1],j[2])
end


## Laplacian 
# Should eventually replace Laplacian and laplacian from Polys.jl
function laplacian(S::ProjectiveMatrixSpace,f)
    sum([derivative(f,gen(S,a,a)) for a in 0:S.N-1])
end

## Ideals

function Ihplus(S::ProjectiveMatrixSpace)
    ideal(S.RP,[hplus(S,j1,j2) for j1 in 0:S.N-1 for j2 in 0:j1])
end

Ihminus(S::ProjectiveMatrixSpace) = ideal(S.RP,[hminus(S,i,j) for i in 0:S.N-1 for j in 0:i if i != j])

function Im(S::ProjectiveMatrixSpace) 
    ideal(S.RP,minors(matrix(S),2))
end

Ireal(S::ProjectiveMatrixSpace) = ideal(S.RP,[gen(S,i,j) - gen(S,j,i) for i in 0:S.N-1 for j in 0:i if i != j])

function It(S::ProjectiveMatrixSpace,t::zzModMatrix) 
    N = S.N
    ideal(S.RP,[gen(S,ZN(N)[j1 j2]) - gen(S,ZN(N)[j1 j2]*t) for j1 in 0:S.N-1 for j2 in 0:S.N-1])
end

function irrelevant_ideal(R::MPolyRing)
    ideal(gens(R))
end

function irrelevant_ideal(S::ProjectiveMatrixSpace)
    irrelevant_ideal(S.RP)
end

function saturation(I::MPolyIdeal)
    saturation(I,irrelevant_ideal(base_ring(I)))
end

function is_saturated(I::MPolyIdeal)
    saturation(I) == I
end