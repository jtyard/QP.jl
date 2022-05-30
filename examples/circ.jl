# Construction of some basic quaternion orders in julia

using Oscar
using QP

# Construct the standard rational quaternion algebra with i^2 = j^2 = k^2 = ijk = -1
A = QuaternionAlgebra(QQ,-1,-1)

# Lets give the standard quaternions names
A1, Ai, Aj, Ak = basis(A)

# Lipschitz quaternion order 
L = Order(A,basis(A)) 
println(discriminant(L))

# Hurwitz quaternion order ∪ ∞
H = Order(A,[A1, Ai, Aj, A(2)^-1*(1 + Ai + Aj + Ak)])
println(discriminant(H))

# Check that H is indeed a maximal order
OA = MaximalOrder(H) # returns a maximal order containing the argument
println(OA == H)


# To build the quaternion order underlying Clifford+T, we need to work over Q(√2)
# To write √, type \sqrt TAB (in vscode with the julia language extension).  
# Horray for unicode ⊗ ⊕ α β γ δ ζ

ZZx, x = PolynomialRing(ZZ, "x")
K, s = NumberField(x^2 - 2, "√2")

B = QuaternionAlgebra(K,-1,-1)
B1, Bi, Bj, Bk = basis(B)

LOK = Order(B,basis(B))
println(discriminant(LOK))

HOK = Order(B,[B1, Bi, Bj, (1 + Bi + Bj + Bk)*K(1//2)])  
println(discriminant(HOK))

# HOK is no longer maximal after extending scalars because
# of my favorite maximal order, which has trivial discriminant

# OBwrong = Order(B,[(1//s)*(B1+Bk), (1//s)*(B1+Bi), (1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2)])
# println(discriminant(OBwrong))
# 
OB = Order(B,[B1, (1//s)*(B1+Bi), (1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2)])
println(discriminant(OB))

# OB == MaximalOrder(OB)
# ρ = α |0⟩⟨0| ⊗ σ + β |1⟩⟨1| ⊗ ω |Ψ⟩

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



    