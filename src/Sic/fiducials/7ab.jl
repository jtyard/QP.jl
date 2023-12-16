using Oscar, QP

# Next do all primes of the form m^2 + 3

S = SicData(7)
K = S.K
OK = S.OK
s = gen(S.K)

# These will work for the order Z[2s] of discriminant 32
# For now (d=7) this gives both fiducials but presumably I'll need these from orders 
# for the others.  Maybe that works already???
# 
#p1 = OK(1+2s)*OK
#p2 = OK(1-2s)*OK

p1 = OK(3-s)*OK
p2 = OK(3+s)*OK

(inf2,inf1) = infinite_places(K) # by convention these are in increasing order so s > 0 in inf1.

rcf1 = ray_class_field(p1,[inf1])
rcf2 = ray_class_field(p2,[inf2])


########
# Scott-Grassl 7a
########
rcfa = cyclotomic_extension(rcf2,4)

Li = number_field(rcfa)

ii = zetaN(4,Li)

e1,e2 = [e for e in complex_embeddings(Li) if real(e(Li(s))) > 0 && imag(e(ii)) > 0]

w = -s
r1 = sqrt(Li(2w+1))

CC = ((-w+2)*r1+(3*w-2))*ii-w*r1-w
AA = ((w+2)*r1+(w+4))*ii+(w-2)*r1-w
BB = 4

A = AA/CC
B = BB/CC

OLi = maximal_order(Li)
a = s*A
b = s*B

psia = Li[1; A; A; B; A; B; B]

c = complex_conjugation(rcfa,inf2)
Psia = psia*dagger(psia,c)


#########
# scott-Grassl 7b
#########

L = number_field(rcf2)

rb = sqrt(L(2s-1))
u = (-rb - 1 - s)/2

A3 = u
B3 = u^-1


psib = L[1;u;u;u^-1;u;u^-1;u^-1]

Psib = psib * transpose(psib)





