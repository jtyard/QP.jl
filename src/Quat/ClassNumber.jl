# Class number and related methods for quaternion orders

using Oscar, QP  

Ord = Hecke.AlgAssRelOrd # could do Union{Hecke.AlgAssAbsOrd,Hecke.AlgAssRelOrd} but disc not an ideal

import Oscar.discriminant, Oscar.prime_divisors, Oscar.class_number

export Ord, split_real_places, ramified places#, prime_divisors, discriminant, class_number

function split_real_places(A::Hecke.QuaternionAlgebra)
    K = base_ring(A)
    a,b = A.std
    return [v for v in real_places(K) if is_negative(a,v) | is_negative(b,v)]
end 

prime_divisors(a::NumFieldElem) = [p[1] for p in factor(a*maximal_order(parent(a))) ]

function ramified_primes(A::Hecke.QuaternionAlgebra)
    OK = maximal_order(base_ring(A))
    a,b = A.std
    [p for p in _candidate_ramified_primes(A) if hilbert_symbol(a,b,p) == -1]
end

function discriminant(A::Hecke.QuaternionAlgebra)
    prod(ramified_primes(A))
end

class_number(A::Hecke.QuaternionAlgebra) = class_number(maximal_order(A))

# Following Section 5 of https://arxiv.org/abs/0808.3833
function _class_number_totally_definite(O::Ord)
    D = discriminant(O)
    facD = factor(D)
end



function class_number(O::Ord)
    A = algebra(O)
    @assert A isa Hecke.QuaternionAlgebra
    K = base_ring(A)
    srp = split_real_places(A)
    if length(srp) == degree(K)
        return _class_number_totally_definite(O)
    else
        return order(ray_class_group(1*maximal_order(K),srp)[1])
    end
end

