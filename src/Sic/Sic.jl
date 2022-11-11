using Oscar
#using Caching


# These are all terrible names to be global in QP.jl
# ultimately will move them out of global after things are working
export minors, Im, QQXp
export XX, h, h2, hp, hm
export Ih, Ih2, Ihm, Ihp, Ic, Icc, Itr0, IT


export SicData


# Making new rings for the overlaps - maybe use routines from Vars.jl instead? 
#Aij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = gens(CA(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
#Aij(j::nmod_mat) = Aij(j[1],j[2],Int(characteristic(base_ring(j))))
#Aij(N::Int) = matrix(CA(N),[[Aij(i,j,N) for j =0:N-1] for i =0:N-1])
#
#aij(j::nmod_mat) = tr(Aij(Int(characteristic(base_ring(j))))*heis(j))










# Not actually used anywhere, subsumed by SICdata
function sicrcf(N::Int)
    D = (N-3)*(N+1)
    _,x = PolynomialRing(QQ)
    K, _ = NumberField(x^2 - D, "s")
    OK = maximal_order(K)
    _,ph = ray_class_group((isodd(N) ? N : 2*N)*OK,infinite_places(K))
    number_field(ray_class_field(ph))
end

# Returns automorphism of F corresponding to place inf of K
function complex_conjugation_from_inf(F,inf)
   automorphism_group(F)
end

mutable struct SicData
    d::fmpz
    D::fmpz    
    D0::fmpz
    f::fmpz
    K::NumField
    OK::NumFieldOrd # the maximal order Z[uf]
    uf::NumFieldOrdElem #or NumFieldElem?
    OD::NumFieldOrd # the minimal order Z[b]
    Ob::NumFieldOrd # the minimal order Z[b]
    b::NumFieldElem # or NumFieldOrdElem?
    rcf::ClassField
    #ring_class_field::ClassField
    function SicData(d::Int;build_nf=false)
        D = ZZ((d-3)*(d+1))
        D0 = fundamental_discriminant(D)
        f = ZZ(sqrt(D//D0))
        K = my_quadratic_field(D0)
        OK = maximal_order(K)
        uf = fundamental_unit(OK)
        OD = quadratic_order(D)
        b = (D-1 + sqrt(K(D)))//2
        Ob = quadratic_order(b)
        rcf = ray_class_field((isodd(d) ? d : 2*d)*OK,infinite_places(K))
        if build_nf
            number_field(rcf)
        end
        #F = number_field(rcf,using_stark_units = true)
        new(d,D,D0,f,K,OK,uf,OD,Ob,b,rcf)
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

Ih(N::Int, graded = false) = ideal(QQX(N, graded = graded),[h(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ih2(N::Int, graded = false) = ideal(QQX(N, graded = graded),[h2(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihp(N::Int, graded = false) = ideal(QQX(N, graded = graded),[hp(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihm(N::Int, graded = false) = ideal(QQX(N, graded = graded),[hm(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])

Im(N::Int, graded = false) = ideal(QQX(N, graded = graded),minors(N))

Ic(N::Int, graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(j2,j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
Icc(N::Int, graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(-j2,-j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
#Itr1(N) = ideal([TrX(N) - 1])
Itr0(N::Int, graded = false) = ideal(QQX(N, graded = graded),[TrX(N)])
IT(a,N::Int, graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2) - Xij(a*j1,a*j2) for j1 in 0:N-1 for j2 in 0:N-1] )

# duh Itr0(N) = R(X_{0,0} + X_{1,1} + ⋯ ) subset RXp(N) = R({X_{i,j}})




