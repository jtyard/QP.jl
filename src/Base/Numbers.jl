# Useful rings, fields and elements
using Oscar

export ZN, zetaN, qint, dagger, primes

#################
# Ring   Element type
# ZZ     fmpz 
# QQ     fmpq
# ZN(N)  nmod
#
# Note that characteristic(base_ring(___)) works for GF(p)[ ] and ResidueClassRing(ZZ,4)
#################

Base.Int(a::nmod) = Int(ZZ(a))
ZN(N) = ResidueRing(ZZ,N) 
zetaN(N) = cyclotomic_field(N)[2]

# Create row vectors over R (or do I want column vectors??)
import Base.^
R^n = MatrixSpace(R,1,n)

# note that rand(R^n) actually works

# conjugate_transpose

dagger(M) = transpose( map(complex_conjugation(base_ring(M)), M ) )

# Perhaps rewrite this to live in the real cyclotomic field
# Returns [m] at q^{1/2} = 2nth root of unity
qint(n,m) = sum([zetaN(2*n)^(m-1 - 2*i) for i in 0:m-1])

function sicrcf(N::Int)
    D = (N-3)*(N+1)
    _,x = PolynomialRing(QQ)
    K, _ = NumberField(x^2 - D, "s")
    OK = maximal_order(K)
    _,ph = ray_class_group((isodd(N) ? N : 2*N)*OK,infinite_places(K))
    number_field(ray_class_field(ph),using_stark_units = true)
end

import Oscar.primes

primes(N::Int) = [n for n in 2:N if is_prime(n)]

primes(M::Int,N::Int) = [n for n in M:N if is_prime(n)]/9