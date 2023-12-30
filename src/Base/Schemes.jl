

using Oscar

import Oscar.gen, Oscar.matrix, Oscar.minors

export ProjectiveMatrixSpace, projective_matrix_space, gen, matrix

export fij, fplus, fminus, hij, hplus, hminus

export Iminors, Ihplus


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

function gen(S::ProjectiveMatrixSpace,i::Union{Int,nmod},j::Union{Int,nmod})
    #@assert divides(Int(S.N),Int(characteristic(base_ring(j))))[1]
    matrix(S)[Int(ZN(S.N)(i))+1,Int(ZN(S.N)(j))+1]
end

function fij(S::ProjectiveMatrixSpace,j1::Union{Int,nmod},j2::Union{Int,nmod}) 
    X(i,j) = gen(S,i,j)
    sum([X(a,a+j1)*X(a+j1+j2,a+j2) for a=0:S.N-1])
end

function fij(S::ProjectiveMatrixSpace,j::nmod_mat) 
    fij(S,j[1],j[2])
end

function fplus(S::ProjectiveMatrixSpace,j1::Union{Int,nmod},j2::Union{Int,nmod}) 
    (fij(S,j1,j2) + fij(S,j2,j1))*S.RP(1//2)
end

function fplus(S::ProjectiveMatrixSpace,j::nmod_mat) 
    fplus(S,j[1],j[2])
end

function fminus(S::ProjectiveMatrixSpace,j1::Union{Int,nmod},j2::Union{Int,nmod}) 
    (fij(S,j1,j2) - fij(S,j2,j1))*S.RP(1//2)
end

function fminus(S::ProjectiveMatrixSpace,j::nmod_mat) 
    fminus(S,j[1],j[2])
end

## Projecting onto harmonic subspace

function hij(S::ProjectiveMatrixSpace,j1::Union{Int,nmod},j2::Union{Int,nmod}) 
    X = gen(S) 
    c = S.RP(1//(2*S.N + 2))*(tr(X^2) + tr(X)^2)
    fij(S,j1,j2) - S.RP((Int(j1) == 0 ? 0 : 1) + (Int(j2) == 0 ? 0 : 1))*c
end

function hij(S::ProjectiveMatrixSpace,j::nmod_mat) 
    hij(S,j[1],j[2])
end

function hplus(S::ProjectiveMatrixSpace,j1::Union{Int,nmod},j2::Union{Int,nmod})
    X = gen(S) 
    c = S.RP(1//(2*S.N + 2))*(tr(X^2) + tr(X)^2)
    0*fplus(S,j1,j2) - S.RP((Int(j1) == 0 ? 0 : 1) + (Int(j2) == 0 ? 0 : 1))*c
end

function hplus(S::ProjectiveMatrixSpace,j::nmod_mat) 
    hplus(S,j[1],j[2])
end

function hminus(S::ProjectiveMatrixSpace,j1::Union{Int,nmod},j2::Union{Int,nmod}) 
    hij(S,j1,j2) - hij(S,j2,j1)
end

function hminus(S::ProjectiveMatrixSpace,j::nmod_mat) 
    hminus(S,j[1],j[2])
end


## Laplacian 
# Should eventually replace Laplacian and laplacian from Polys.jl
function laplacian(S::ProjectiveMatrixSpace,f)
    sum([derivative(f,gen(S,a,a)) for a in 0:S.N-1])
end

## Ideals

function Ihplus(S)
    ideal(S.RP,[hplus(S,j1,j2) for j1 in 0:S.N-1 for j2 in 0:j1])
end


function Im(S) 
    ideal(S.RP,minors(matrix(S),2))
end