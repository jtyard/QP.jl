using Oscar

include("Sic.jl")

N = 7

#CN = cyclotomic_field(N)

ZZx, x = PolynomialRing(ZZ, "x")
D = (N-3)*(N+1)

K, s = NumberField(x^2 - 2, "s")
OK = maximal_order(K)

infK = infinite_places(K)

if real(evaluate(s,infK[1])) < 0
    infK = [infK[2],infK[1]]
end

PP = prime_decomposition(OK,N)
P1 = PP[1][1]
P2 = PP[2][1]

G1, ph1 = ray_class_group(P1,[infK[2]])
G2, ph2 = ray_class_group(P2,[infK[2]])

n1 = order(G1)
n2 = order(G2)
if n2 < n1
    P1,G1,ph1,n1,P2,G2,ph2,n2 = P2,G2,ph2,n2,P1,G1,ph1,n1
end

Ky, y = PolynomialRing(K,"y")
F,r = NumberField(y^2 - (2*s-1),"r")
# Make it work for ii later

ZZz,zs = PolynomialRing(ZZ,[string("z",i) for i in 0:N-1])
z(i) = zs[1 + (i % N)]

norm2 = sum([z(i)^2 for i in 0:N-1])
norm4 = norm2^2

zzzz(j) = sum([z(i)*z(i+j[1])*z(i+j[2])*z(i+j[1]+j[2]) for i in 0:N-1])
f(j) = (N+1)*zzzz(j) - norm4*( (((j[1] % N) == 0) ? 1 : 0) + (((j[2] % N) == 0) ? 1 : 0))

u = (r-1-s)//2 
psi = [F(1),u,u,u^-1,u,u^-1,u^-1]

[[evaluate(f([j1,j2]),psi) for j2 in 0:N-1] for j1 in 0:N-1]

CCC() = AcbField(precision(BigFloat))
RRR() = ArbField(precision(BigFloat))

function starks(N::Integer)
    L = lvals(N)
    n = length(L)
    CC = CCC()
    Larb = [CC(real(a)) + onei(CC)*CC(imag(a)) for a in L]
    zn = exppii(CC(2/n))
    zetas = [sum([zn^(i*j) * Larb[j] for j = 1:2:n]) for i = 0:n-1]
    [exp(real(a)) for a in zetas]
end



function findquadpoly(a::arb)
    N = ZZ(2)^200
    MS = MatrixSpace(ZZ,3,4)
    A = MS([1 0 0 N; 0 1 0 ZZ(floor(N*a)); 0 0 1 ZZ(floor(N*a^2))])
    println(A)
    B = lll(A)
    B[1,1] + B[1,2]*x+B[1,3]*x^2
end

function arbstarkcoefs(N::Integer)
    _,t = PolynomialRing(RRR(),"t")
    f = prod([t-s for s in starks(N)])
    [a for a in coefficients(f)]
end

# Quick and dirty weil reprsentation for diagonals
function wA(a::Integer,N::Integer)
    U = zeros(Integer,N,d)
    U[1,1] = 1
    for j = 1:N-1
        U[j + 1,(a*j % N) + 1] = 1

        
    end
    #U*jacobi_symbol(a,d)
    U # ignore the jacobi sign to always fix |0>
end




