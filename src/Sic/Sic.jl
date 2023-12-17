# Definitions relevant to studying SICs.  SicData(d) constructs a number of objects relevant to SICs in dimension d.  SicData(d,build_nf=true) will also construct the ray class field. 
# Harmonic invariant polynomials labeled by elements of (Z/d)^2


using Oscar
using Memoize
#using Caching

include("SicData.jl")
include("Harm2.jl")

# These are all terrible names to be global in QP.jl
# ultimately will move them out of global after things are working
export minors, Im, QQXp
export XX, h, h2, hp, hm
export Ih, Ih2, Im, Ihm, Ihp, Ic, Icc, Itr0, Itr1, IT

import Oscar.overlaps, Oscar.complex_conjugation
export overlaps, complex_conjugation
export is_fiducial, heis_orbit, algebra_from_heis_orbit, algebra_from_basis

export dwork_modulus




function is_fiducial(Phi::AbstractAlgebra.Generic.MatSpaceElem) 
    N = ncols(Phi)
    F = base_ring(Phi)
    for j in 0:N-1
        for i in 0:N-1
            if evaluate(h(ZN(N)[i j]),vec(Phi)) != F(0)
                return false
            end
        end
    end
    return true
end
                
function overlaps(Phi::AbstractAlgebra.Generic.MatSpaceElem) 
    N = ncols(Phi)
    F = base_ring(Phi)
    C_to_F = hom(cyclotomic_field(N)[1],F,zetaN(N,F))
    matrix(F,N,N,[trace(map(C_to_F,heis(ZN(N)[i j]))*Phi) for j in 0:N-1 for i in 0:N-1])
end


function heis_orbit(Phi)
    N = ncols(Phi)
    NN = (isodd(N) ? N : 2*N)
    F = base_ring(Phi)
    C_to_F = hom(cyclotomic_field(NN)[1],F,zetaN(NN,F))
    [map(C_to_F,heis(ZN(N)[i j]))*Phi* map(C_to_F,heis(ZN(N)[i j])^-1) for j in 0:N-1 for i in 0:N-1]
end


function sic(d::Int)
    heis_orbit(fiducial(d))
end

# e.g. fiducial("7b") returns the Scott-Grassl 7b sic
function sic(label::String)
    heis_orbit(fiducial(label))
end

function algebra_from_basis(Phis)
    F = base_ring(Phis[1])
    matrix_algebra(QQ,F,Phis)
end


function algebra_from_heis_orbit(Phi)
    F = base_ring(Phi)
    A = matrix_algebra(QQ,F,heis_orbit(Phi))
    AA = AlgAss(A)[1]
    AZ = Hecke._as_algebra_over_center(AA)[1]
end



function is_split(A::AlgMat)
    is_split(AlgAss(A)[1])
end


# SICs on Calabi-Yaus
# Verifying that each fiducial lies on a unique Dwork hypersurface X_1^N + ... + X_N^N - N mu X_1 ... X_N
# with e^{t/N} a unit in the complex ring class field of the corresponding quadratic order.
# 
# In particular this should be real when there are real fiducials.

function dwork_modulus(P)
    if nrows(P) == 1 || ncols(P) == 1
        v = P
    else
        v = P[:,1]
    end
    N = length(v)
    sum([a^N for a in v])//prod([a for a in v]) # Often divided by N but it is more integral without doing that
end