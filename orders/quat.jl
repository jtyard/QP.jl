# Construction of some basic quaternion orders in julia

using Oscar
using QP


include("su2k.jl")

# Lipschitz quaternion order 
function lipschitz_quaternions()
    # Construct the standard rational quaternion algebra with i^2 = j^2 = k^2 = ijk = -1
    A = quaternion_algebra(QQ,-1,-1)
    A1, Ai, Aj, Ak = basis(A)
    Order(A,basis(A)) 
end


# Hurwitz quaternion order
function hurwitz_quaternions()
    A = quaternion_algebra(QQ,-1,-1)
    A1, Ai, Aj, Ak = basis(A)
    H = Order(A,[A1, Ai, Aj, A(2)^-1*(1 + Ai + Aj + Ak)])
end

# Check that H is indeed a maximal order

println(hurwitz_quaternions() == MaximalOrder(hurwitz_quaternions())

# To build the quaternion order underlying Clifford+T, we need to work over Q(√2)
# "Clifford" quaternion order
function clifford_quaternions()
    
    ZZx, x = PolynomialRing(ZZ, "x")
    K, s = NumberField(x^2 - 2, "√2") # \sqrt TAB = \sqrt 


    B = quaternion_algebra(K,-1,-1)
    B1, Bi, Bj, Bk = basis(B)

    #write this into lipschitz_quaternions(OK)
    #LOK = Order(B,basis(B))
    #println(discriminant(LOK))

    #HOK = Order(B,[B1, Bi, Bj, (1 + Bi + Bj + Bk)*K(1//2)])  
    #println(discriminant(HOK))

    # HOK is no longer maximal after extending scalars because
    # of my favorite maximal order, which has trivial discriminant

    # OBwrong = Order(B,[(1//s)*(B1+Bk), (1//s)*(B1+Bi), (1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2)])
    # println(discriminant(OBwrong))
    # 

    OB = Order(B,[B1, (1//s)*(B1+Bi), (1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2)])
end
    println(discriminant(OB))

    # OB == MaximalOrder(OB)

bb = [B1, (1//s)*(B1+Bi),(1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2), s*B1, B1+Bi, B1+Bj, (B1+Bi+Bj+Bk)*(s//2)]
bbconj =[B1, (1//s)*(B1-Bi),(1//s)*(B1-Bj), (B1-Bi-Bj-Bk)*K(1//2), s*B1, B1-Bi, B1-Bj, (B1-Bi-Bj-Bk)*(s//2)]
bbconjprime =[B1, -(1//s)*(B1-Bi),-(1//s)*(B1-Bj), (B1-Bi-Bj-Bk)*K(1//2), -s*B1, B1-Bi, B1-Bj, -(B1-Bi-Bj-Bk)*(s//2)]
bbprime =[B1, -(1//s)*(B1+Bi),-(1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2), -s*B1, B1+Bi, B1+Bj, -(B1+Bi+Bj+Bk)*(s//2)]

G1 = matrix([[trace(trace(a*b)//(8+4*s)) for b in bb] for a in bb])
G2 = matrix([[trace(trace(a*b)//(8+4*s)) for b in bbconj] for a in bb])



println(det(G1))
println(det(G2))
println(eigvals(G1))
println(eigvals(G2))


su2 = su2levelk(3)


# Need numerics

RRR() = ArbField(precision(BigFloat))
CCC() = AcbField(precision(BigFloat))
ii = onei(CCC())

RR = RRR()
CC(x) = CCC()(x) 
CC(M::Array) = map(x->CC(x),M)

eigvals(matrix([[CC(trace(trace(a*b)//(8+4*s))) for b in bb] for a in bbconj]))



# assuming splitting at K(i)
function toM(q)
    K = q.parent.base_ring
    _,x = PolynomialRing(K)
    N = normred(q)
    F,i = NumberField(x^2 +1,"i")
    P = [matrix(F,[1 0; 0 1]), matrix(F,[0 1; 1 0]), matrix(F,[0 -i; i 0]), matrix(F,[1 0; 0 -1])]

    P[1]*F(q.coeffs[1]) - i * sum([P[j]*F(q.coeffs[j]) for j in 2:4]) 
end

function toU(q) 
    M = toM(q)
    F = M.base_ring
    N = det(M)

    if issquare(N)[1]
        sqrtN = sqrt(N)
        MM = M
    else
        _,x = PolynomialRing(F)
        FF,sqrtN = NumberField(x^2 - N)
        MM = matrix(FF,M)
    end
    MM
end



    
