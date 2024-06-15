# Exact synthesis and catalytic embeddings for qutrits

using Oscar


Q3, z3 = cyclotomic_field(3)
s3 = z3 - z3^-1

Q9 = cyclotomic_extension(Q3,9).Kr

z9 = gen(Q9)
chi = z9 - z9^-1

H = Q9[1 1 1; 1 z3 z3^2; 1 z3^2 z3]/s3

S = Q9[1 0 0; 0 1 0; 0 0 z3]

T = Q9[1 0 0; 0 z9 0; 0 0 z9^-1]

R = Q9[1 0 0; 0 1 0; 0 0 -1] # "metaplectic gate" outside Clifford hierarchy


