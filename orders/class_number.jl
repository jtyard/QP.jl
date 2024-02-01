# Class number and related methods for quaternion orders

using Oscar, QP  


K,s = quadratic_field(2)
A = quaternion_algebra(QQ,-1,-1)
B = quaternion_algebra(K,1,1+s)

OA = maximal_order(A)
OB = maximal_order(B)


Ord = Hecke.AlgAssRelOrd # could do Union{Hecke.AlgAssAbsOrd,Hecke.AlgAssRelOrd} but disc not an ideal



class_number(A::Hecke.QuaternionAlgebra) = class_number(maximal_order(A))
split_real_places(A::Hecke.QuaternionAlgebra) = split_real_places(maximal_order(A))


# Following Section 5 of https://arxiv.org/abs/0808.3833
function _class_number_totally_definite(O::Ord)
    D = discriminant(O)
    facD = factor(D)
end

function split_real_places(A::Hecke.QuaternionAlgebra)
    K = base_ring(A)
    a,b = A.std
    return [v for v in real_places(K) if is_negative(a,v) | is_negative(b,v)]
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

