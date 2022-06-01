#################
# Algebras and orders
#################
export QuaternionAlgebra, Dn, su2levelk

# A wrapper that coerces the generators automatically 
QuaternionAlgebra(K,a,b) = Hecke.AlgQuat(K,K(a),K(b))

# From paper with Vadym for SU(2)_k
export Dn, su2levelk
Dn(n) = n % 4 == 0 ? 1 : 4*(1-CyclotomicRealSubfield(n)[2]^2)

# q = zeta_{k+2}, q^{1/2} = zeta_{2k + 4}

# The quaternion algebra (K, -[3]_q,-D_n) 
function su2levelk(k) 
    K, qq = CyclotomicRealSubfield(k+2)
    QuaternionAlgebra(K,-1-qq,-Dn(k+2)) 
end