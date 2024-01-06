# Useful rings, fields and elements
using Oscar, Memoize

export Z2, Z3, Z4, Z5, Z6, Z7, Z8, ZN

export zetaN, qint, dagger, complex_conjugate, primes, myorder

export roots_of_unity, primitive_roots_of_unity

#################
# Ring   Element type
# ZZ     ZZRingElem
# QQ     QQFieldElem
# ZN(N)  zzModRingElem 
#
# Note that characteristic(base_ring(___)) works for GF(p)[ ] and ResidueClassRing(ZZ,4)
#################

Base.Int(a::zzModRingElem) = Int(ZZ(a))

# Caching works unexpectedly - `@memoize` ensures that e.g. Z2 == ZN(2) always, which was broken for `using QP`
@memoize ZN(N) = ResidueRing(ZZ,N) 
Z2 = ZN(2)
Z3 = ZN(3)
Z4 = ZN(4)
Z5 = ZN(5)
Z6 = ZN(6)
Z7 = ZN(7)
Z8 = ZN(8)
zetaN(N) = cyclotomic_field(N)[2]

roots_of_unity(F,n) = roots(polynomial_ring(F)[2]^n-1)
primitive_roots_of_unity(F,n) = roots(polynomial_ring(F)[1](collect(coefficients(cyclotomic_field(n)[1].pol))))

zetaN(n,F) = primitive_roots_of_unity(F,n)[1]

# Maybe rewrite this to live in the real cyclotomic field
# Returns [m] at q^{1/2} = 2nth root of unity
function qint(n,m) 
    if m==0 
        return 0
    end
    q = iseven(n) ? zetaN(2*n)  : -zetaN(n)
    sum([q^(m-1 - 2*i) for i in 0:m-1])
end

# Another approach (used to be in Fmatrix.jl)
#_, q = RationalFunctionField(QQ,"q")
#qint(n) = sum([q^(n-1 - 2*i) for i in 0:n-1])


# Create row vectors over R (or do I want column vectors??)
import Base.^
R^n = MatrixSpace(R,n,1)

# note that rand(R^n) actually works

# conjugate_transpose

#dagger(M::MatElem,c) = transpose( map(x->map(c,x), M))
dagger(M::MatElem,c) = transpose(map(c,M))
dagger(M::MatElem) = dagger(M,complex_conjugation(base_ring(M)))

complex_conjugate(x) = complex_conjugation(parent(x))(x)


import Oscar.primes

primes(N::Int) = collect(PrimesSet(2,N))

primes(M::Int,N::Int) = collect(PrimesSet(M,N))

# quick and dirty in case order(a) doesn't work 
function myorder(a)
    aa = a
    id = a^0
    for i in 1:1000
        if aa == id
            return i
        end
        aa = aa*a
    end
    return Inf
end 

