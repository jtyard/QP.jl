# Playing with gram matrices

using Oscar, QP

OB = clifford_quaternions()
B = algebra(OB)
B1, Bi, Bj, Bk = basis(B)
K = base_ring(B)
s = gen(K)

bb = [B1, (1//s)*(B1+Bi),(1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2), s*B1, B1+Bi, B1+Bj, (B1+Bi+Bj+Bk)*(s//2)]
bbconj =[B1, (1//s)*(B1-Bi),(1//s)*(B1-Bj), (B1-Bi-Bj-Bk)*K(1//2), s*B1, B1-Bi, B1-Bj, (B1-Bi-Bj-Bk)*(s//2)]
bbconjprime =[B1, -(1//s)*(B1-Bi),-(1//s)*(B1-Bj), (B1-Bi-Bj-Bk)*K(1//2), -s*B1, B1-Bi, B1-Bj, -(B1-Bi-Bj-Bk)*(s//2)]
bbprime =[B1, -(1//s)*(B1+Bi),-(1//s)*(B1+Bj), (B1+Bi+Bj+Bk)*K(1//2), -s*B1, B1+Bi, B1+Bj, -(B1+Bi+Bj+Bk)*(s//2)]

G1 = matrix([[trace(trace(a*b)//(8+4*s)) for b in bb] for a in bb])
G2 = matrix([[trace(trace(a*b)//(8+4*s)) for b in bbconj] for a in bb])

println(det(G1))
println(det(G2))
println(eigvals(G1))
println(eigvals(G2))

# Floating point representation

RRR() = ArbField(precision(BigFloat))
CCC() = AcbField(precision(BigFloat))
ii = onei(CCC())

RR = RRR()
CC(x) = CCC()(x) 
CC(M::Array) = map(x->CC(x),M)

eigvals(matrix([[CC(trace(trace(a*b)//(8+4*s))) for b in bb] for a in bbconj]))