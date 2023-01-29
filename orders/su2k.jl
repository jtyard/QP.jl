#################
# Algebras and orders for SU(2)_k
#################
export Dn, su2levelk, su2levelkD, su2levelkR

# Some orders from https://arxiv.org/abs/1504.04350 (Section 4) and https://arxiv.org/abs/1510.03888 for SU(2)_k

# Is totally positive and sqrt(-Dn(n)) generates Q(zeta_n) 
Dn(n) = n % 4 == 0 ? 1 : 4-CyclotomicRealSubfield(n)[2]^2

# q = zeta_{k+2}, q^{1/2} = zeta_{2k + 4}

# The quaternion algebra (K, -[3]_q,-D_n) 
function su2levelk(k) 
    K, qq = CyclotomicRealSubfield(k+2)
    quaternion_algebra(K,-1-qq,-Dn(k+2)) 
end

su2levelkD(k) = discriminant(maximal_order(su2levelk(k)))

su2levelkR(k) = base_ring(maximal_order(su2levelk(k)))