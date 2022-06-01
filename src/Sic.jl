using Oscar
using Caching


#################
# Ring   Element type
# ZZ     fmpz 
# QQ     fmpq
# ZN(N)  nmod
#
# Can I use this yet??? https://github.com/oscar-system/Oscar.jl/blob/master/test/Schemes/ProjectiveSchemes.jl
#
# Note that characteristic(base_ring(___)) works for GF(p)[ ] and ResidueClassRing(ZZ,4)
#################

Base.Int(a::nmod) = Int(ZZ(a))

ZN(N) = ResidueRing(ZZ,N) 

#M(N) = FreeModule(ZN(N),2)
#gl(N) = GL(2,ZN(N))
#sl(N) = SL(2,ZN(N))

#################
# Matrix units and projectors
#################

# Matrix units
Eij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = matrix(ZZ,[[a == (Int(i) % N) && b == (Int(j) % N) for b in 0:N-1] for a in 0:N-1])

# Access with an nmod vector like Eij(ZN(N)[i j])
Eij(ij::nmod_mat) = Eij(ij[1], ij[2], Int(characteristic(base_ring(ij))))

# 2d array of matrix units 
Eij(N::Int) = [[Eij(i,j,N) for j in 0:N-1] for i in 0:N-1]



SWAP(N) = sum([kronecker_product(Eij(i,j,N),Eij(j,i,N)) for i in 1:N for j in 1:N])
Pplus(N) = (1//2)*(identity_matrix(QQ,N^2) + SWAP(N))
Pminus(N) = (1//2)*(identity_matrix(QQ,N^2) - SWAP(N))


###################
# Matrix polynomials
###################

# Oscar only caches rings not graded rings - I need QQX(N) to always refer to the *same* ring
@cache function QQX(N::Int)
    GradedPolynomialRing(QQ,[string("X_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1])[1]
end

Xij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = gens(QQX(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
Xij(j::nmod_mat) = Xij(j[1],j[2],Int(characteristic(base_ring(j))))
Xij(N::Int) = matrix(QQX(N),[[Xij(i,j,N) for j =0:N-1] for i =0:N-1])

@cache function CA(N::Int)
    GradedPolynomialRing(cyclotomic_field(N)[1],[string("A_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1])[1]
end

# Not sure what Oscar can do with these at the moment (i.e. with Proj) but I can ask later. 
@cache function QQzw(N::Int)
    GradedPolynomialRing(QQ,vcat([string("z",i) for i in 0:N-1],[string("w",i) for i in 0:N-1]))[1]
end

zj(j::nmod) = gens(QQX(Int(characteristic(base_ring(j)))))[1 + (i % N)]
wj(j::nmod) = gens(QQX(Int(characteristic(base_ring(j)))))[1 + N + (i % N)]

# Overlaps
Aij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = gens(CA(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
Aij(j::nmod_mat) = Aij(j[1],j[2],Int(characteristic(base_ring(j))))
Aij(N::Int) = matrix(CA(N),[[Aij(i,j,N) for j =0:N-1] for i =0:N-1])

aij(j::nmod_mat) = tr(Aij(Int(characteristic(base_ring(j))))*heis(j))

# Ring homomorphism taking Xij(N) -> Y when Y is an N x N matrix
QQXhom(Y) = hom(QQX(nrows(Y)),QQX(nrows(Y)),vec(transpose(Y)))
QQXt(N::Int) = QQXhom(transpose(Xij(N)))


TrX(N::Int) = sum([Xij(i,i,N) for i = 0:N-1])
TrX2(N::Int) = sum([Xij(i,j,N)*Xij(j,i,N) for i = 0:N-1 for j = 0:N-1])

# The 2x2 minors cut out the rank < 2 matrices
minors(N::Int) = [Xij(i,k,N)*Xij(j,l,N)-Xij(i,l,N)*Xij(j,k,N) for i in 0:N-2 for j in i+1:N-1 for k in 0:N-2 for l in k+1:N-1];

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

# The irrelevant ideal generated by the matrix units
QQXp(N::Int) = ideal(QQX(N),gens(QQX(N)))

# duh Itr0(N) = R(X_{0,0} + X_{1,1} + ⋯ ) subset RXp(N) = R({X_{i,j}})

############
# Number fields
############

function sicrcf(N::Int)
    D = (N-3)*(N+1)
    _,x = PolynomialRing(QQ)
    K, _ = NumberField(x^2 - D, "s")
    OK = maximal_order(K)
    _,ph = ray_class_group((isodd(N) ? N : 2*N)*OK,infinite_places(K))
    number_field(ray_class_field(ph),using_stark_units = true)
end



