

using Oscar

import Oscar.gen, Oscar.matrix, Oscar.minors

export ProjectiveMatrixSpace, projective_matrix_space, gen, minors, matrix

export Im, fij, hij



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

function gen(S::ProjectiveMatrixSpace,i::Union{Int,nmod},j::Union{Int,nmod})
    #@assert divides(Int(S.N),Int(characteristic(base_ring(j))))[1]
    matrix(S)[Int(ZN(S.N)(i))+1,Int(ZN(S.N)(j))+1]
end

function fij(S::ProjectiveMatrixSpace,j1,j2) 
    X(i,j) = gen(S,i,j)
    sum([X(a,a+j1)*X(a+j1+j2,a+j2) for a=0:S.N-1])
end

function hplus(S::ProjectiveMatrixSpace,j1,j2) 
    X = matrix(S)
    
end

function Im(S) 
    ideal(S.RP,minors(matrix(S)))
end