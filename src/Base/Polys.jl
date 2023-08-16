###############
# Polynomial rings and their generators in some convenient forms for computing with polynomial functions 
# on matrices and projective spaces
#
# In the future it may be better to start with a ProjectiveSpace or ProjectiveScheme etc and inherit 
# polynomials / graded polynomials from there.
###############

using Oscar
using Memoize 

export MatrixPolynomialRing, VariableMatrix, MatrixGradedPolynomialRing, QQXgraded
export FX, QQX, Xij, TrX, TrX2, QQXhom, QQXt
export QQzw, wj, zj, monomials_of_degree, laplacian, Laplacian
export QQXtozw

###############
# Matrix polynomials 
###############

# Oscar doesn't cache this yet
@memoize function MatrixGradedPolynomialRing(F::AbstractAlgebra.Ring,N::Int;X::Union{AbstractString, Char, Symbol} = "X")
    GradedPolynomialRing(F,[string(X,"_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1],)[1] 
end

function MatrixPolynomialRing(F,N::Int; X::Union{AbstractString, Char, Symbol} = "X", graded = false)
    if graded == false
        PolynomialRing(F,[string(X,"_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1],cached=true)[1]
    else
        MatrixGradedPolynomialRing(F,N,X)
    end
end

function VariableMatrix(F,N::Int;X::Union{AbstractString, Char, Symbol} = "X", graded = false)
    R = MatrixPolynomialRing(F,N,X,graded=graded)
    matrix(R,N,N,gens(R))
end



#####
function FX(F::AbstractAlgebra.Ring,N::Int; graded = false) 
    graded ? MatrixGradedPolynomialRing(F,N) : MatrixPolynomialRing(F,N)
end

function QQX(N::Int; graded = false) 
   FX(QQ,N,graded = graded) 
end

# Okay for now but a more general definition allowing other names and instances ultimately needed
Xij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int; graded = false) = gens(QQX(N,graded=graded))[1 + (Int(j) % N) + N*(Int(i) % N)]
Xij(j::nmod_mat; graded = false) = Xij(j[1],j[2],Int(characteristic(base_ring(j))),graded=graded)
Xij(N::Int; graded = false) = matrix(QQX(N,graded=graded),[[Xij(i,j,N,graded=graded) for j =0:N-1] for i =0:N-1])


# Ring homomorphism taking Xij(N) -> Y when Y is an N x N matrix
QQXhom(Y) = hom(QQX(nrows(Y)),QQX(nrows(Y)),vec(transpose(Y)))
QQXt(N::Int) = QQXhom(transpose(Xij(N)))

# Mapping to z and w 
QQXtozw(N) = hom(QQX(N),QQzw(N),[gens(QQzw(N))[i]*gens(QQzw(N))[N+j] for i=1:N for j=1:N])


TrX(N::Int; graded = false) = sum([Xij(i,i,N,graded=graded) for i = 0:N-1])
TrX2(N::Int; graded = false) = sum([Xij(i,j,N,graded=graded)*Xij(j,i,N,graded=graded) for i = 0:N-1 for j = 0:N-1])

# The 2x2 minors cut out the rank < 2 matrices
minors(N::Int; graded = false) = [Xij(i,k,N,graded=graded)*Xij(j,l,N,graded=graded)-Xij(i,l,N,graded=graded)*Xij(j,k,N,graded=graded) for i in 0:N-2 for j in i+1:N-1 for k in 0:N-2 for l in k+1:N-1];

# The irrelevant ideal generated by the matrix units
QQXp(N::Int; graded = false) = ideal(QQX(N,graded=graded),gens(QQX(N,graded=graded)))

# Accepts a matrix in SL(2,Z/N) or SL(2,Z/2N) and acts on QQX(N)
#function weil_X(g::nmod_mat)
#    N = Int(characteristic(base_ring(g)))
#    if divides(N,4)[1]
        



#####################
# z zbar coordinates 
#####################

@memoize function QQzw(N::Int; graded = false)
    graded ? GradedPolynomialRing(QQ,vcat([string("z",i) for i in 0:N-1],[string("w",i) for i in 0:N-1]), vcat([ [1,0] for i in 1:N],[[0,1] for i in 1:N]))[1] : PolynomialRing(QQ,vcat([string("z",i) for i in 0:N-1],[string("w",i) for i in 0:N-1]))[1]
end

# Assuming for now input is of the type QQzw - can/should be generalized
function monomials_of_degree(R::MPolyRing_dec,n::Union{Int,fmpz})
    Rnn, to_R = homogeneous_component(R,[n,n])
    [to_R(f) for f in gens(Rnn)]
end

function zj(j::nmod) 
    N = Int(characteristic(parent(j)))
    j = Int(j)
    gens(QQzw(N))[1 + (j % N)]
end

function wj(j::nmod) 
    N = Int(characteristic(parent(j)))
    j = Int(j)
    gens(QQzw(N))[1 + N + (j % N)]
end






laplacian(f) = sum([derivative(derivative(f,zj(a)),wj(a)) for a in ZN(ZZ(length(gens(parent(f)))//2))])

function Laplacian(f)
    N = Int(sqrt(length(gens(parent(gens(Ih(4))[1])))))
    sum([derivative(f,Xij(a,a,N)) for a in 0:N-1])
end

# Not sure Oscar is being consistent with their usage of "base_ring" as they should sometimes call it "ambient ring"
import Oscar.change_base_ring

function change_base_ring(F::AbstractAlgebra.Ring,R::MPolyRing_dec)
    GradedPolynomialRing(F,[string(x) for x in R.R.S],R.d)[1]
end

function change_base_ring(F::AbstractAlgebra.Ring,R::MPolyRing)
    PolynomialRing(F,[string(x) for x in R.S])
end

function change_base_ring(F::AbstractAlgebra.Ring,I::MPolyIdeal)
    R = base_ring(I)
    RF = change_base_ring(F,R)
    ideal(RF,[change_base_ring(F,f,parent=RF) for f in gens(I)])
end


