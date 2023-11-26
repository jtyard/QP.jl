# Definitions relevant to studying SICs.  SicData(d) constructs a number of objects relevant to SICs in dimension d.  SicData(d,build_nf=true) will also construct the ray class field. 
# Harmonic invariant polynomials labeled by elements of (Z/d)^2


using Oscar
using Memoize
#using Caching


export SicData

# These are all terrible names to be global in QP.jl
# ultimately will move them out of global after things are working
export minors, Im, QQXp
export XX, h, h2, hp, hm
export Ih, Ih2, Im, Ihm, Ihp, Ic, Icc, Itr0, Itr1, IT

import Oscar.overlaps, Oscar.complex_conjugation
export overlaps, complex_conjugation
export is_fiducial, heis_orbit, algebra_from_heis_orbit, algebra_from_basis, fiducial, sic


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
    F::NumField
    c::Hecke.NumFieldMor
    e::Hecke.NumFieldEmbNfNS
    #ring_class_field::ClassField
    @memoize function SicData(d::Int;build_nf=true)
        if d == 3
            return new(3,0,0,1,rationals_as_number_field()[1])
        end
        D = ZZ((d-3)*(d+1))
        D0 = fundamental_discriminant(D) 
        f = ZZ(sqrt(D//D0)) 
        K = quadratic_field(Hecke.squarefree_part(D0))[1] 
        inf = real_places(K)
        OK = maximal_order(K)
        uf = d > 3 ? fundamental_unit(OK) : OK(1)
        OD = quadratic_order(D)
        bb = (d-1 + sqrt(K(D)))//K(2)
        Ob = quadratic_order(bb)
        b = Ob(bb)
        rcf = ray_class_field((isodd(d) ? d : 2*d)*OK,infinite_places(K))
    
        S = new(d,
        D,
        D0,
        f,
        K,
        inf,
        OK,
        uf,
        OD,
        Ob,
        b,
        rcf)
        if build_nf
            S.F = number_field(rcf)
            F = S.F
            #S.F = number_field(rcf,using_stark_units = true)
            S.c = complex_conjugation(rcf,inf[2])
            
            # would be nice if we could do
            # complex_conjugation(S.F) = c 
            # but it's not so easily possible

            CC = AcbField(64) #bitsize should be increased for larger fiducials but this will work for small ones. Revisit this.
            zdd = exp(2*pi*onei(CC)/(2*d))
            S.e = [e for e in complex_embeddings(F) if overlaps(e(zetaN(2*d,F)),zdd) && overlaps(e(sqrt(F(D0))),sqrt(CC(D0)))][1]
        end
        S
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
    matrix(F,N,N,[trace(map(C_to_F,heis(ZN(N)[i j]))*Phi) for j in 0:N-1 for i in 0:N-1])
end


function heis_orbit(Phi)
    N = ncols(Phi)
    NN = (isodd(N) ? N : 2*N)
    F = base_ring(Phi)
    C_to_F = hom(cyclotomic_field(NN)[1],F,zetaN(NN,F))
    [map(C_to_F,heis(ZN(N)[i j]))*Phi* map(C_to_F,heis(ZN(N)[i j])^-1) for j in 0:N-1 for i in 0:N-1]
end

function algebra_from_heis_orbit(Phi)
    F = base_ring(Phi)
    A = matrix_algebra(QQ,F,heis_orbit(Phi))
    AA = AlgAss(A)[1]
    AZ = Hecke._as_algebra_over_center(AA)[1]
end

# Returns a fiducial in dimension d - todo: return more like SICdata, rcf, complex conjugation?
function fiducial(d::Int)
    if d == 2
        F = cyclotomic_field(12)[1]
        return (identity_matrix(F,2) + map(F,heis(0,1,2) + heis(1,0,2) + heis(1,1,2))/sqrt(F(3)))/2
    elseif d == 3
        F = cyclotomic_field(3)[1]
        return F[0;1;-1]*F[0 1 -1]/2
    elseif d == 4
        S4 = SicData(4,build_nf=true)
        F = S4.F
        #F = number_field(S4.rcf,using_stark_units = true)

        # Define some constants
        w2 = sqrt(F(2))
        w5 = sqrt(F(5))
        w10 = w2*w5
        r1=sqrt(w5+1)
        I=sqrt(F(-1))
        phi=S4.F[
            8
            ((w10+w2-2*w5-2)*r1-4)*I+(w10+w2)*r1+4*w2-4
            (8*w2-8)*I
        ((w10+w2-2*w5-2)*r1+4)*I+(w10+w2)*r1-4*w2+4]

        #gal, sig = automorphism_group(S4.rcf)
        #abs2c(x) = x*c(x)
     
        # Rank-1 density matrix 
        Phi = phi*dagger(phi,S4.c); 
        return trace(Phi)^-1 * Phi
    elseif d == 5
        S5 = SicData(5, build_nf=true)
        F = S5.F

        # Define some constants
        w3 = sqrt(F(3));
        w5 = sqrt(F(5));
        w15 = w3*w5;
        I = sqrt(F(-1));

        r1=sqrt(1//2*w5+5//2)
        r2=sqrt((-1//8*w3+1//8*w5+3//8)*r1+1//8*w3*w5+5//8*w3+1//16*w5-5//16)

        # Find a complex conjugation 
        #gal, sig = automorphism_group(S5.rcf.A)
        #possible_c = findall([map(sig(g),[I,w3,w5,r1,r2]) == [-I,w3,w5,r1,r2] for g in gal])
        #println("Better only be one choice in this list " * string(possible_c))
        #c = sig(gal[possible_c[1]])
        #gal, sig = automorphism_group(S5.rcf)
        #c = complex_conjugation(S5.rcf, S5.inf[2])
        #abs2c(x) = x*c(x)

        # The fiducial vector
        phi = F[16*r1
            (((8*w3-10*w5-2)*r1+20*w3-16*w5)*r2+((2*w5+2)*r1+(2*w5+10)))*I+((-6*w15-6*w3+8*w5+24)*r1+(-12*w15-20*w3+12*w5+40))*r2+(4*w5-8)*r1-4*w15-4*w5
            (((16*w3-14*w5+6)*r1+(-2*w15+30*w3-24*w5))*r2+((w15-5*w3+4)*r1-3*w5+5))*I+((2*w15-2*w3+4*w5-12)*r1-2*w5-10)*r2+(w5+3)*r1+3*w15-5*w3-2*w5+10
            (((-6*w15-6*w3+6*w5+14)*r1+(-8*w15-20*w3+16*w5+20))*r2+((w15+5*w3-2*w5+6)*r1+4*w15-2*w5))*I+((2*w15+2*w3-2*w5-18)*r1+20*w3-12*w5)*r2+(3*w5-1)*r1-2*w15
            (((-8*w15+12*w3-12*w5+32)*r1+(-10*w15+10*w3-6*w5+30))*r2+((-2*w15+w5+3)*r1+w15+5*w3+3*w5+5))*I+((-4*w3-4)*r1+2*w15-10*w3+2*w5-10)*r2+(-w15+5*w3+2*w5)*r1+3*w15+5*w3+w5+5] 

        # Rank-1 density matrix 
        Phi = phi*dagger(phi,S5.c); 
        return trace(Phi)^-1 * Phi
    else
        error("Not implemented")
    end
end

# e.g. fiducial("7b") returns the Scott-Grassl 7b fiducial
function fiducial(label::String)
    if label == "7a"
        error("Not implemented")
    elseif label == "7b"
        error("Not implemented")
    else
        error("Not implemented")
    end
end

function sic(d::Int)
    heis_orbit(fiducial(d))
end

# e.g. fiducial("7b") returns the Scott-Grassl 7b sic
function sic(label::String)
    error("Not implemented")
end

function algebra_from_basis(Phis)
    F = base_ring(Phis[1])
    matrix_algebra(QQ,F,Phis)
end


function is_split(A::AlgMat)
    is_split(AlgAss(A)[1])
end