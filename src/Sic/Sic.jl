using Oscar
#using Caching


# These are all terrible names to be global in QP.jl
# ultimately will move them out of global after things are working
export minors, Im, QQXp
export XX, h, h2, hp, hm
export Ih, Ih2, Ihm, Ihp, Ic, Icc, Itr0, IT

export my_real_embeddings, my_quadratic_field, quad_embedding, fundamental_unit, quadratic_order

export SICdata, SICdata_nf


# Making new rings for the overlaps - maybe use routines from Vars.jl instead? 
#Aij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = gens(CA(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
#Aij(j::nmod_mat) = Aij(j[1],j[2],Int(characteristic(base_ring(j))))
#Aij(N::Int) = matrix(CA(N),[[Aij(i,j,N) for j =0:N-1] for i =0:N-1])
#
#aij(j::nmod_mat) = tr(Aij(Int(characteristic(base_ring(j))))*heis(j))







#######################
# Quadratic number fields  
#########################

# Real embeddings should be actually real
# cf https://github.com/thofma/Hecke.jl/issues/491
my_real_embeddings(K::NumField) = [real ∘ e for e in real_embeddings(K)]

# real embeddings of quadratic_field(m) have tiny imaginary parts.  Doing this instead.
function my_quadratic_field(m::Union{fmpz,Int})
    _,x = PolynomialRing(QQ)
    NumberField(x^2-m,'s')[1]
end

# The quadratic order of discriminant D
function quadratic_order(D::Union{fmpz,Int})
    @assert mod(D,4) in [0,1] string(D) * " ≢ 0,1 mod 4}"
    K = my_quadratic_field(fundamental_discriminant(D))
    NfAbsOrd([K(1), (K(D) + sqrt(K(D)))//2])
end

# The quadratic order of discriminant D
function quadratic_order(b::NumFieldElem)
    K = parent(b)
    @assert degree(K) == 2 "Ambient field of " * string(b) * " must be quadratic"
    @assert degree(b) == 2 string(b) * " must be quadratic"
    NfAbsOrd([K(1), b])
end

# Choose the embdding making the generator positive 
function quad_embedding(K::NumField) 
    embs = my_real_embeddings(K)
    embs[indexin(1,[e(gen(K)) > 0 for e in embs])[1]]
end

# Fundamental unit is smallest unit > 1 
function fundamental_unit(OK::NumFieldOrd) 
    u = unit_group(OK)[2](unit_group(OK)[1]([0,1]))
    K = OK.nf
    e = quad_embedding(K)
    if e(K(u)) < 0
        u = -u
    end
    if e(K(u)) < 1
        u = u^-1
    end
    return u
end

# Fundamental unit of the maximal order
fundamental_unit(K::NumField) = fundamental_unit(maximal_order(K))

# Fundamental unit of the quadratic order.
fundament_unit(D::Union{fmpz,Int}) = fundamental_unit(quadratic_order(D))


# Not actually used anywhere, subsumed by SICdata
function sicrcf(N::Int)
    D = (N-3)*(N+1)
    _,x = PolynomialRing(QQ)
    K, _ = NumberField(x^2 - D, "s")
    OK = maximal_order(K)
    _,ph = ray_class_group((isodd(N) ? N : 2*N)*OK,infinite_places(K))
    number_field(ray_class_field(ph),using_stark_units = true)
end

# Returns automorphism of F corresponding to place inf of K
function complex_conjugation_from_inf(F,inf)
   automorphism_group(F)
end

mutable struct SICdata
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
    function SICdata(d::Int)
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
        #F = number_field(rcf,using_stark_units = true)
        new(d,D,D0,f,K,OK,uf,OD,Ob,b,rcf)
    end
end

function SICdata_nf(d::Int)
    S = SICdata(d)
    number_field(S.rcf)
    return S
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

Ih(N::Int) = ideal(QQX(N),[h(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ih2(N::Int) = ideal(QQX(N),[h2(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihp(N::Int) = ideal(QQX(N),[hp(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihm(N::Int) = ideal(QQX(N),[hm(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])

Im(N::Int) = ideal(QQX(N),minors(N))

Ic(N::Int) = ideal(QQX(N), [Xij(j1,j2,N) - Xij(j2,j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
Icc(N::Int) = ideal(QQX(N), [Xij(j1,j2,N) - Xij(-j2,-j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
#Itr1(N) = ideal([TrX(N) - 1])
Itr0(N::Int) = ideal(QQX(N),[TrX(N)])
IT(a,N::Int) = ideal(QQX(N), [Xij(j1,j2) - Xij(a*j1,a*j2) for j1 in 0:N-1 for j2 in 0:N-1] )

# duh Itr0(N) = R(X_{0,0} + X_{1,1} + ⋯ ) subset RXp(N) = R({X_{i,j}})




