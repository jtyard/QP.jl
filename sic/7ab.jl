using Oscar, QP

# Get this working for 7 then for all primes of the form m^2 + 3


S = SicData(7)
K = S.K
OK = S.OK
#_,y = polynomial_ring(K)
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



#LL,(r,rr,ii) = number_field([y^2 - (2*s -1),y^2 - (-2*s-1),y^2+1])

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


art = artin_map(rcfa)

G = domain(art.map2)
gal(g) = art.map2(G(g))

println( gal([1,1])(a) == b)
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

is_fiducial(Psia)

#LL<r,rr,ii> := ext<K|y^2 - (2*s -1),y^2 - (-2*s-1),y^2+1>;
#s7 := r*rr;
#r1 := ii*rr;

#F := RelativeField(K,LL);
#G := Automorphisms(F);

#aut := func<k| map<LL->LL |x:-> (LL ! G[k](F ! x))> >;
#temp := [[Integers() ! (aut(k)(x)/x) : x in [r,rr,ii]]:k in [1..8]];


#c1 := aut(Index(temp,[1,-1,-1]));
#c2 := aut(Index(temp,[-1,1,-1]));
#sig1 := aut(Index(temp,[-1,1,1]));
#sig2 := aut(Index(temp,[1,-1,1]));


# scott-Grassl 7a:
#CC := ((-s+2)*r1+(3*s-2))*ii-s*r1-s;
#AA := ((s+2)*r1+(s+4))*ii+(s-2)*r1-s;
#    BB := 4;

#z8 := (1+ii)/s;

#A := AA/CC;
#B := BB/CC;        

println(discriminant(maximal_order(cyclotomic_field(8)[1])))

mu = dwork_modulus(psia)*14
is_isomorphic(cyclotomic_field(8)[1],number_field(absolute_minpoly(mu))[1])

