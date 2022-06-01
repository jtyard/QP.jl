###############
# Polynomial rings and their generators in some convenient forms for computing with polynomial functions 
# on matrices and projective spaces
###############

export MatrixPolynomialRing, VariableMatrix
export QQX,Xij,QQzw,wj,zj

###############
# Polynomials 
###############

# Generators 

function MatrixPolynomialRing(F,N::Int,X::Union{AbstractString, Char, Symbol} = "X")
    PolynomialRing(F,[string(X,"_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1],cached=true)[1]
end

function VariableMatrix(F,N::Int,X::Union{AbstractString, Char, Symbol} = "X")
    R = MatrixPolynomialRing(F,N,X)
    matrix(R,N,N,gens(R))
end

#####
# "Legacy"
#####
QQX(N) = MatrixPolynomialRing(QQ,N)

Xij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = gens(QQX(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
Xij(j::nmod_mat) = Xij(j[1],j[2],Int(characteristic(base_ring(j))))
Xij(N::Int) = matrix(QQX(N),[[Xij(i,j,N) for j =0:N-1] for i =0:N-1])


# Not sure what Oscar can do with these at the moment (i.e. with Proj) but I can ask later. 
@cache function QQzw(N::Int)
    GradedPolynomialRing(QQ,vcat([string("z",i) for i in 0:N-1],[string("w",i) for i in 0:N-1]))[1]
end

zj(j::nmod) = gens(QQX(Int(characteristic(base_ring(j)))))[1 + (i % N)]
wj(j::nmod) = gens(QQX(Int(characteristic(base_ring(j)))))[1 + N + (i % N)]