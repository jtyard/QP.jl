using Oscar, QP


# Returns the affine cone in C^{d^n} over the projective scheme of perfect tensors of type (d,n)
# If t is set it restricts to the real sphere S^{2d^n-1} of squared radius t
function perfect_tensors(d::Int, n::Int; t = 0)
    R = MatPolyRing(QQ,d^n)
    X = gen(R);

    k = floor(Int,n/2)
    subs = subsets(n,k)
    rhos = [reduced_operator(X,[d for _ in 1:n],s) for s in subs]

    ideals =[ideal(R,vec(collect(d^k*rho - trace(rho)*diagonal_matrix(R(1),d^k)))) for rho in rhos]

    S = polynomial_ring(QQ,vcat([Symbol("z_", i) for i in 0:d^n-1],[Symbol("zbar_",j) for j in 0:d^n-1]))[1];
    _phi = hom(R,S,[gen(S,i)*gen(S,j+d^n) for i in 1:d^n for j in 1:d^n]);

    I = _phi(sum(ideals) + ideal(R,t != 0 ? tr(X)-R(t) : R(0)))

    spec(I)
end

# Constructs a rational point of the variety with coordinates zi = zbari = p[i].
# Assumes the p[i] are real as it's all we need so far - can add complex conjugation if needed
# Needs better error handling as ideals are so huge the error leaves the buffer
# In such case for now can use  write("error_log.txt", sprint(showerror,err)) to see the last error
function _point(P,p) 
    p = [p[i] for i in 1:length(p)] # need to flatten p for ket(..) to work as it's a matrix.
    #c = sum([a^2 for a in p])
    P(vcat(p,p))  
end

cat_state(d,n) = sum([ket(matrix(ZN(d),transpose([i for _ in 1:n]))) for i in  0:d-1])

# computes dimension of tangent space at the point p 
dim_at_point(P,p) = dim(tangent_space(_point(P,p)))


println(dim_at_point(perfect_tensors(2,3), cat_state(2,3)))
# 9 
println(dim_at_point(perfect_tensors(2,3,t=2), cat_state(2,3)))
# 8
# projective would cut out the U(1) and give 7


# Rational points 
# d = 3, n = 3
per = ket(Z3[0 1 2]) + ket(Z3[1 2 0]) + ket(Z3[2 0 1]) + ket(Z3[2 1 0]) + ket(Z3[1 0 2]) + ket(Z3[0 2 1]);
ghz3 = ket(Z3[0 0 0]) + ket(Z3[1 1 1]) + ket(Z3[2 2 2]);
_det= ket(Z3[0 1 2]) + ket(Z3[1 2 0]) + ket(Z3[2 0 1]) - ket(Z3[2 1 0]) - ket(Z3[1 0 2]) - ket(Z3[0 2 1]);

P = perfect_tensors(3,3)
println((dim_at_point(P,per), dim_at_point(P,cat_state(3,3)), dim_at_point(P,_det)))
# (32, 34, 38)


# Could try (2,2,2,2,2), (3,3,3,3) or (4,4,4) but I need a rational point
#@time P3333 = perfect_tensors(3,4)
#p = _point(P44,___) # says it's not a rational point
#@time T = tangent_space(p)
#@time dim(T)
