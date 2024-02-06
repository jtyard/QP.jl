# Class number and related methods for quaternion orders

using Oscar, QP  

#import Hecke.discriminant, Oscar.prime_divisors, Hecke.QuaternionAlgebra, Oscar.NumFieldElem

export Ord

export split_real_places, class_number, ramified_places, class_number, prime_divisors, discriminant, class_number

Ord = Hecke.AlgAssRelOrd # could do Union{Hecke.AlgAssAbsOrd,Hecke.AlgAssRelOrd} but disc not an ideal
#QuaternionAlgebra  = Hecke.QuaternionAlgebra

# A wrapper that coerces the generators automatically - should merge into OSCAR
quaternion_algebra(K,a,b) = Hecke.QuaternionAlgebra(K,K(a),K(b))

function split_real_places(A::Hecke.QuaternionAlgebra)
    K = base_ring(A)
    a,b = A.std
    return [v for v in real_places(K) if is_negative(a,v) | is_negative(b,v)]
end 

Oscar.prime_divisors(a::NumFieldElem) = [p[1] for p in factor(a*maximal_order(parent(a))) ]

function Oscar.ramified_primes(A::Hecke.QuaternionAlgebra)
    K = base_ring(A)
    a,b = A.std
    [p for p in union(prime_divisors(a),prime_divisors(b),prime_divisors(K(2))) if hilbert_symbol(a,b,p) == -1]
end

function Oscar.discriminant(A::Hecke.QuaternionAlgebra)
    OK = maximal_order(base_ring(A))
    prod([ramified_primes(A);1*OK])
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

