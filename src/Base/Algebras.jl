#################
# Algebras and orders
#################
export quaternion_algebra, Dn, su2levelk, su2levelkD, su2levelkR

# A wrapper that coerces the generators automatically 
quaternion_algebra(K,a,b) = Hecke.AlgQuat(K,K(a),K(b))

# From paper with Vadym for SU(2)_k

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