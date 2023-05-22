#################
# Algebras and orders for SU(2)_k
#################
export Dn, su2levelk

# Some orders from https://arxiv.org/abs/1504.04350 (Section 4) and https://arxiv.org/abs/1510.03888 for SU(2)_k

# Dn(n) is totally positive and sqrt(-Dn(n)) generates Q(zeta_n):
Dn(n) = n % 4 == 0 ? 1 : 4-CyclotomicRealSubfield(n)[2]^2

# q = zeta_{k+2}, q^{1/2} = zeta_{2k + 4}

mutable struct su2levelk
    K::AnticNumberField
    OK::NfOrd
    qq::nf_elem # 2nth root of unity
    A::Hecke.AlgQuat # The quaternion algebra (K, -[3]_q,-D_n)
    OA::Hecke.AlgAssRelOrd 
    D::NfOrdIdl
 
    function su2levelk(k::Int) 
        K, qq = CyclotomicRealSubfield(k+2)
        A = quaternion_algebra(K,-1-qq,-Dn(k+2)) 
        OK = maximal_order(K)
        OA = maximal_order(A)
        D = discriminant(OA)
        new(K,OK,qq,A,OA,D)
    end
end

#su2levelkD(k) = discriminant(maximal_order(su2levelk(k)))

#su2levelkR(k) = base_ring(maximal_order(su2levelk(k)))