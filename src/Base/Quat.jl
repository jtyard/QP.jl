# Construction of some basic quaternion orders in julia

using Oscar
using QP

export quaternion_algebra, lipschitz_quaternions, hurwitz_quaternions, clifford_quaternions, su2levelk, toM, toU

# A wrapper that coerces the generators automatically - should merge into OSCAR
quaternion_algebra(K,a,b) = Hecke.AlgQuat(K,K(a),K(b))

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
function clifford_quaternions()
    ZZx, x = PolynomialRing(ZZ, "x")
    K, s = NumberField(x^2 - 2, "√2") # \sqrt TAB = \sqrt 
    clifford_quaternions(K)
end

function clifford_quaternions(K::Field)
    @assert is_square_with_sqrt(K(2))[1]
    s = sqrt(K(2))

    B = quaternion_algebra(K,-1,-1)
    B1, Bi, Bj, Bk = basis(B)

    Order(B,[B1, (1//s)*(B1+Bi), (1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2)])
end

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

    
#println(discriminant(OB))

# OB == MaximalOrder(OB)


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

# Some orders from https://arxiv.org/abs/1504.04350 (Section 4) and https://arxiv.org/abs/1510.03888 for SU(2)_k

# Dn(n) is totally positive and sqrt(-Dn(n)) generates Q(zeta_n):
#Dn(n) = n % 4 == 0 ? 1 : 4-CyclotomicRealSubfield(n)[2]^2

# q = zeta_{k+2}, q^{1/2} = zeta_{2k + 4}

mutable struct su2levelk
    K::AnticNumberField
    qq::nf_elem # 2nth root of unity
    Dn::nf_elem # Dn(n) is totally positive and sqrt(-Dn(n)) generates Q(zeta_n):
    OK::NfOrd
    A::Hecke.AlgQuat # The quaternion algebra (K, -[3]_q,-D_n)
    OA::Hecke.AlgAssRelOrd 
    D::NfOrdIdl
 
    function su2levelk(k::Int) 
        n = k+2
        K, qq = CyclotomicRealSubfield(n)
        Dn = n % 4 == 0 ? 1 : 4-CyclotomicRealSubfield(n)[2]^2
        A = quaternion_algebra(K,-1-qq,-Dn) 
        OK = maximal_order(K)
        OA = maximal_order(A)
        D = discriminant(OA)
        new(K,qq,Dn,OK,A,OA,D)
    end
end

    
