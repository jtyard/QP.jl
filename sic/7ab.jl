using Oscar, QP

# Get this working for 7 then for all primes of the form m^2 + 3


S = SicData(7)
K = S.K
OK = S.OK
s = gen(S.K)


LL<r,rr,ii> := ext<K|y^2 - (2*s -1),y^2 - (-2*s-1),y^2+1>;
s7 := r*rr;
r1 := ii*rr;

F := RelativeField(K,LL);
G := Automorphisms(F);

aut := func<k| map<LL->LL |x:-> (LL ! G[k](F ! x))> >;
temp := [[Integers() ! (aut(k)(x)/x) : x in [r,rr,ii]]:k in [1..8]];


c1 := aut(Index(temp,[1,-1,-1]));
c2 := aut(Index(temp,[-1,1,-1]));
sig1 := aut(Index(temp,[-1,1,1]));
sig2 := aut(Index(temp,[1,-1,1]));


#_<z> := PolynomialRing(L);
#L2<ii,s7> := ext<L|z^2 + 1,z^2 + 7>;


# scott-Grassl 7a:
CC := ((-s+2)*r1+(3*s-2))*ii-s*r1-s;
AA := ((s+2)*r1+(s+4))*ii+(s-2)*r1-s;
    BB := 4;

z8 := (1+ii)/s;

A := AA/CC;
B := BB/CC;        