# Definitions relevant to studying SICs.  SicData(d) constructs a number of objects relevant to SICs in dimension d.  SicData(d,build_nf=true) will also construct the ray class field. 
# Harmonic invariant polynomials labeled by elements of (Z/d)^2


using Oscar
#using Caching


export SicData

# These are all terrible names to be global in QP.jl
# ultimately will move them out of global after things are working
export minors, Im, QQXp
export XX, h, h2, hp, hm
export Ih, Ih2, Im, Ihm, Ihp, Ic, Icc, Itr0, Itr1, IT

import Oscar.overlaps
export is_fiducial, overlaps, heis_orbit, algebra_from_basis


# Making new rings for the overlaps - maybe use routines from Vars.jl instead? 
#Aij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = gens(CA(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
#Aij(j::nmod_mat) = Aij(j[1],j[2],Int(characteristic(base_ring(j))))
#Aij(N::Int) = matrix(CA(N),[[Aij(i,j,N) for j =0:N-1] for i =0:N-1])
#
#aij(j::nmod_mat) = tr(Aij(Int(characteristic(base_ring(j))))*heis(j))

# Not actually used anywhere, subsumed by SICdata

mutable struct SicData
    d::fmpz
    D::fmpz    
    D0::fmpz
    f::fmpz
    K::NumField
    inf::Vector{InfPlc}
    OK::NumFieldOrd # the maximal order Z[uf]
    uf::NumFieldOrdElem #or NumFieldElem?
    OD::NumFieldOrd # the minimal order Z[b]
    Ob::NumFieldOrd # the minimal order Z[b]
    b::NumFieldOrdElem # or NumFieldOrdElem?
    rcf::ClassField
    #ring_class_field::ClassField
    function SicData(d::Int;build_nf=false)
        D = ZZ((d-3)*(d+1))
        D0 = fundamental_discriminant(D)
        f = ZZ(sqrt(D//D0))
        K = quadratic_field(Hecke.squarefree_part(D0))[1]
        inf = real_places(K)
        OK = maximal_order(K)
        uf = fundamental_unit(OK)
        OD = quadratic_order(D)
        bb = (D-1 + sqrt(K(D)))//2
        Ob = quadratic_order(bb)
        b = Ob(bb)
        rcf = ray_class_field((isodd(d) ? d : 2*d)*OK,infinite_places(K))
        if build_nf
            number_field(rcf)
        end
        #F = number_field(rcf,using_stark_units = true)
        new(d,D,D0,f,K,inf,OK,uf,OD,Ob,b,rcf)
    end
end

X = Xij

function XX(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    sum([X(a,a+j[1],N)*X(a+j[1]+j[2],a+j[2],N) for a in 0:N-1])
end

XX(j1,j2,N::Int) = XX(ZN(N)[j1 j2])

function h(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    (N+1)*XX(j) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*TrX(N)^2
end
h(j1,j2,N::Int) = h(ZN(N)[j1 j2])

function h2(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    (N+1)*XX(j) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*TrX2(N)
end
h2(j1,j2,N::Int) = h2(ZN(N)[j1 j2])

function hp(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    ((N+1)//2)*(XX(j) + XX(ZN(N)[j[2] j[1]])) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*(TrX(N)^2 + TrX2(N))*(1//2)
end
hp(j1,j2,N::Int) = hp(ZN(N)[j1 j2])

function hm(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    ((N+1)//2)*(XX(j) - XX(ZN(N)[j[2] j[1]])) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*(TrX(N)^2 - TrX2(N))*(1//2)
end
hm(j1,j2,N::Int) = hm(ZN(N)[j1 j2])


############
# Useful ideals
############

Ih(N::Int; graded = false) = ideal(QQX(N, graded = graded),[h(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ih2(N::Int; graded = false) = ideal(QQX(N, graded = graded),[h2(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihp(N::Int; graded = false) = ideal(QQX(N, graded = graded),[hp(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihm(N::Int; graded = false) = ideal(QQX(N, graded = graded),[hm(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])

Im(N::Int; graded = false) = ideal(QQX(N, graded = graded),minors(N))

Ic(N::Int; graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(j2,j1,N) for j1 in 0:N-1 for j2 in 0:N-1 if j1 < j2] )
Icc(N::Int; graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(-j2,-j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
Itr1(N) = ideal([TrX(N) - 1])
Itr0(N::Int; graded = false) = ideal(QQX(N, graded = graded),[TrX(N)])
IT(a,N::Int; graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(ZN(N)(a)*j1,ZN(N)(a)*j2,N) for j1 in 0:N-1 for j2 in 0:N-1] )

# duh Itr0(N) = R(X_{0,0} + X_{1,1} + ⋯ ) subset RXp(N) = R({X_{i,j}})



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
    [[trace(map(C_to_F,heis(ZN(N)[i j]))*Phi) for j in 0:N-1] for i in 0:N-1]
end


function heis_orbit(Phi)
    N = ncols(Phi)
    NN = (isodd(N) ? N : 2*N)
    F = base_ring(Phi)
    C_to_F = hom(cyclotomic_field(NN)[1],F,zetaN(NN,F))
    [map(C_to_F,heis(ZN(N)[i j]))*Phi* map(C_to_F,heis(ZN(N)[i j])^-1) for j in 0:N-1 for i in 0:N-1]
end

function algebra_from_basis(Phis)
    F = base_ring(Phis[1])
    matrix_algebra(QQ,F,Phis)
end


function is_split(A::AlgMat)
    is_split(AlgAss(A)[1])
end