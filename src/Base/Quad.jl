using Oscar

export quaternion_algebra, my_real_embeddings, my_quadratic_field, quad_embedding, fundamental_unit, quadratic_order

#######################
# Useful things for quadratic number fields, orders, etc
# Add some to Oscar?
#########################

# A wrapper that coerces the generators automatically - should merge into OSCAR
quaternion_algebra(K,a,b) = Hecke.AlgQuat(K,K(a),K(b))

# Real embeddings should be actually real
# cf https://github.com/thofma/Hecke.jl/issues/491
my_real_embeddings(K::NumField) = [real ∘ e for e in real_embeddings(K)]

# real embeddings of quadratic_field(m) have tiny imaginary parts.  Doing this instead.
function my_quadratic_field(m::Union{fmpz,Int})
    _,x = PolynomialRing(QQ)
    NumberField(x^2-m,'s')[1]
end

# The quadratic order of discriminant D
function quadratic_order(D::Union{fmpz,Int})
    @assert mod(D,4) in [0,1] string(D) * " ≢ 0,1 mod 4}"
    K = my_quadratic_field(fundamental_discriminant(D))
    NfAbsOrd([K(1), (K(D) + sqrt(K(D)))//2])
end

# The quadratic order of discriminant D
function quadratic_order(b::NumFieldElem)
    K = parent(b)
    @assert degree(K) == 2 "Ambient field of " * string(b) * " must be quadratic"
    @assert degree(b) == 2 string(b) * " must be quadratic"
    NfAbsOrd([K(1), b])
end

# Choose the embdding making the generator positive 
function quad_embedding(K::NumField) 
    embs = my_real_embeddings(K)
    embs[indexin(1,[e(gen(K)) > 0 for e in embs])[1]]
end

# Fundamental unit is smallest unit > 1 
function fundamental_unit(OK::NumFieldOrd) 
    u = unit_group(OK)[2](unit_group(OK)[1]([0,1]))
    K = OK.nf
    e = quad_embedding(K)
    if e(K(u)) < 0
        u = -u
    end
    if e(K(u)) < 1
        u = u^-1
    end
    return u
end

# Fundamental unit of the maximal order
fundamental_unit(K::NumField) = fundamental_unit(maximal_order(K))

# Fundamental unit of the quadratic order.
fundamental_unit(D::Union{fmpz,Int}) = fundamental_unit(quadratic_order(D))