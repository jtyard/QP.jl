

# export print_with_subscripts, display_with_subscripts, latex_representation

struct MatPolyRing{T, R<:Ring, S<:MPolyRing} <: Ring
    R::R
    N::Int
    S::S
    Sgr::MPolyDecRing

    function MatPolyRing(R::Ring, N::Int, vars::Vector{String}; internal_ordering=:lex)
        S, _ = polynomial_ring(R, vars; internal_ordering=internal_ordering)
        T = Oscar.elem_type(S)
        new{T, typeof(R), typeof(S)}(R, N, S, grade(S)[1])
    end
end

# Constructor
function MatPolyRing(F::Ring, N::Int; X::Union{AbstractString, Char, Symbol} = "X", internal_ordering=:lex)
    vars = [string(X, "_{", i, ",", j, "}") for i in 0:N-1 for j in 0:N-1]
    MatPolyRing(F, N, vars; internal_ordering=internal_ordering)
end

# Alternate constructor that accepts a symbol for the variable name
#function MatPolyRing(F::Ring, N::Int; X::Symbol = :X)
#    MatPolyRing(F, N, X=string(X))
#end


matrix_polynomial_ring(F::Ring, N::Int; X::Union{AbstractString, Char, Symbol} = "X", internal_ordering=:lex) = MatPolyRing(F, N; X=X, internal_ordering=internal_ordering) 
#matrix_polynomial_ring(F::Ring, N::Int; X::Symbol, internal_ordering=:lex) = MatPolyRing(F, N; X=X, internal_ordering=internal_ordering)





## Display methods

function Base.show(io::IO, R::MatPolyRing)
    X = split(string(R.S.S[1]),'_')[1]
    print(io, "Polynomial ring in $(R.N)×$(R.N) variables $(X)_{⋅,⋅} over ", R.R)
end


# Helper function to convert subscript numbers
function to_subscript(n::Int)
    subscript_digits = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
    return join(subscript_digits[digit(n, i) + 1] for i in reverse(1:ndigits(n)))
end



# Implement necessary methods
Oscar.base_ring(R::MatPolyRing) = R.R
Oscar.ngens(R::MatPolyRing) = R.N^2
Oscar.gens(R::MatPolyRing) = gens(R.S)
Oscar.gen(R::MatPolyRing, i::Int) = gen(R.S, i)

# Add matrix-specific methods
function gen(R::MatPolyRing)
    matrix(R.S, R.N, R.N, gens(R.S))
end

function Oscar.gen(R::MatPolyRing, i::Int, j::Int)
    gen(R.S, (i % R.N)*R.N + (j % R.N) + 1)
end

# Element type
Oscar.elem_type(::Type{<:MatPolyRing{T}}) where {T} = T

# Define zero and one for the ring
function Oscar.zero(R::MatPolyRing)
    return matrix(R.S, R.N, R.N, [zero(R.S) for _ in 1:R.N^2])
end

function Oscar.one(R::MatPolyRing)
    return matrix(R.S, R.N, R.N, [i == j ? one(R.S) : zero(R.S) for i in 1:R.N for j in 1:R.N])
end

# Conversion method
function (R::MatPolyRing)(a)
    return R.S(a)
end

# homomorphisms
function Oscar.hom(R::MatPolyRing, S::Ring, imgs::Vector{<:RingElem})

    return Oscar.MPolyAnyMap(R.S, S, nothing, copy(imgs)) # copy because of #655
end

# Ideal constructors
function Oscar.ideal(R::MatPolyRing, g::Union{RingElem, RingElement, Integer}...)
    ideal(R.S, g...)
end

function Oscar.ideal(R::MatPolyRing, g::Vector)
    ideal(R.S, g)
end

# Quotient ring constructor
function Oscar.quo(R::MatPolyRing, I::MPolyIdeal)
    Q, f = quo(R.S, I)
    return Q, f
end




