# Construction of some basic quaternion orders in julia

using Oscar
using QP

export Ord

export lipschitz_quaternions, hurwitz_quaternions, clifford_quaternions, su2levelk, toM, toU



Ord = Hecke.AlgAssRelOrd # could do Union{Hecke.AlgAssAbsOrd,Hecke.AlgAssRelOrd} but disc not an ideal
#QuaternionAlgebra  = Hecke.QuaternionAlgebra


#include("su2k.jl")

# Lipschitz quaternion order 
function lipschitz_quaternions(F::Field)
    # Construct the standard rational quaternion algebra with i^2 = j^2 = k^2 = ijk = -1
    A = quaternion_algebra(F,-1,-1)
    A1, Ai, Aj, Ak = basis(A)
    Order(A,basis(A)) 
end

lipschitz_quaternions() = lipschitz_quaternions(QQ)

# Hurwitz quaternion order
function hurwitz_quaternions(F::Field)
    A = quaternion_algebra(F,-1,-1)
    A1, Ai, Aj, Ak = basis(A)
    H = Order(A,[A1, Ai, Aj, A(2)^-1*(1 + Ai + Aj + Ak)])
end

hurwitz_quaternions() = hurwitz_quaternions(QQ)

# Check that H is indeed a maximal order
# println(hurwitz_quaternions() == MaximalOrder(hurwitz_quaternions()))

# To build the quaternion order underlying Clifford+T, we need to work over Q(√2)
# "Clifford" quaternion order

function clifford_quaternions(K::Field)
    if is_square(K(2))
        F = K
        s = sqrt(K(2))
    else
        _, x = polynomial_ring(K, "x")
        F, s = NumberField(x^2 - 2, "√2") # \sqrt TAB = \sqrt 
    end

    B = quaternion_algebra(F,-1,-1)
    B1, Bi, Bj, Bk = basis(B)

    Order(B,[B1, (1//s)*(B1+Bi), (1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2)])
end

 

    # OBwrong = Order(B,[(1//s)*(B1+Bk), (1//s)*(B1+Bi), (1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2)])
    # println(discriminant(OBwrong))


clifford_quaternions() = clifford_quaternions(QQ)
    
#println(discriminant(OB))

# OB == MaximalOrder(OB)


# assuming splitting at K(i)
function toM(q)
    K = q.parent.base_ring
    _,x = polynomial_ring(K)
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
        _,x = polynomial_ring(F)
        FF,sqrtN = NumberField(x^2 - N)
        MM = matrix(FF,M)
    end
    MM
end

# Some orders from https://arxiv.org/abs/1504.04350 (Section 4) and https://arxiv.org/abs/1510.03888 for SU(2)_k

# Dn(n) is totally positive and sqrt(-Dn(n)) generates Q(zeta_n):
#Dn(n) = n % 4 == 0 ? 1 : 4-CyclotomicRealSubfield(n)[2]^2

# q = zeta_{k+2}, q^{1/2} = zeta_{2k + 4}

mutable struct su2levelk
    K::AnticNumberField
    qq::NumFieldElem # 2nth root of unity
    Dn::NumFieldElem # Dn(n) is totally positive and sqrt(-Dn(n)) generates Q(zeta_n):
    OK::NumFieldOrder
    A::Hecke.QuaternionAlgebra # The quaternion algebra (K, -[3]_q,-D_n)
    OA::Hecke.AlgAssRelOrd 
    D::NumFieldOrderIdeal
 
    function su2levelk(k::Int) 
        n = k+2
        K, qq = CyclotomicRealSubfield(n)
        Dn = n % 4 == 0 ? K(1) : 4-qq^2
        A = quaternion_algebra(K,-1-qq,-Dn) 
        OK = maximal_order(K)
        OA = maximal_order(A)
        D = discriminant(OA)
        new(K,qq,Dn,OK,A,OA,D)
    end
end

    
