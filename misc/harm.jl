@time using Oscar

N = 5;

RX,Xvars = Singular.PolynomialRing(QQ,[string("X_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1])
X(j) = Xvars[1 + (j[2] % N) + N*(j[1] % N)]


R,vars = Singular.PolynomialRing(QQ,vcat([string("z",i) for i in 0:N-1],[string("w",i) for i in 0:N-1]))
z(i) = vars[1 + (i % N)]
w(i) = vars[1 + N + (i % N)]
z(j::Array) = prod([z(i-1)^j[i] for i in 1:N])
w(j::Array) = prod([w(i-1)^j[i] for i in 1:N])

norm2 = sum([z(i)*w(i) for i in 0:N-1])
norm4 = norm2^2

zwwz(j) = sum([z(i)*w(i+j[1])*w(i+j[2])*z(i+j[1]+j[2]) for i in 0:N-1])
f(j) = (N+1)*zwwz(j) - norm4*( (((j[1] % N) == 0) ? 1 : 0) + (((j[2] % N) == 0) ? 1 : 0))
f(j1,j2) = f([j1,j2])

J = [[j1,j2] for j1 in 0:N-1 for j2 in 0:j1]

carprod(A,B) = [vcat(a,b) for a in A for b in B]
function carpow(A,n)
    B = A 
    for i in 1:n-1
        B = carprod(A,B)
    end
    B
end

indofdegree(N,n) = [j for j in carpow(0:n,N) if sum(j) == n]

monomialsofdegree(n) = [z(j)*w(k) for j in indofdegree(N,n) for k in indofdegree(N,n) ];

ZN = ResidueRing(ZZ,N)

function coeffofmonomial(f,m)
    inds = findall(x->x==m, collect(monomials(f)))
    if length(inds) == 0
        0
    else
        collect(coeffs(f))[inds[1]]
    end
end

function tovec(f,n)
    # Assuming f is homogeneous of degree (n,n)
    M = monomialsofdegree(n)
    [QQ(coeffofmonomial(f,m)) for m in M]
end

tovec(f) = tovec(f, ZZ(total_degree(f)/2))

M2 = monomialsofdegree(2)
M3 = monomialsofdegree(3)

Rn(n) = VectorSpace(QQ,length(monomialsofdegree(n)))
R2 = Rn(2)

R3 = Rn(3)

H,phiH = sub(R2,[R2(tovec(f(j))) for j in J])
dim(H)

Del(f) = sum([derivative(derivative(f,z(i)),w(i)) for i in 0:N-1])
Deln(n) = hcat([tovec(Del(m),n-1) for m in monomialsofdegree(n)]...)
Hn(n) = sub(Rn(n),[Rn(n)(tovec(Del(m),n-1) for m in monomialsofdegree(n))])
Deln(2)