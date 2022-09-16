using Oscar
using QP


Dn(n) = n % 4 == 0 ? 1 : 4*(1-CyclotomicRealSubfield(n)[2]^2)

# From paper with Vadym for SU(2)_k
# q = zeta_{k+2}, q^{1/2} = zeta_{2k + 4}
function su2levelk(k) 
    K, qq = CyclotomicRealSubfield(k+2)
    QuaternionAlgebra(K,-1-qq,-Dn(k+2)) # (K, -[3]_q,-D_n)
end
