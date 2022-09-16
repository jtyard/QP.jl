// Build the Scott-Grassl 7a fiducial
// Note that it already has the symmetry I want

d:=7;

_<x> := PolynomialRing(Rationals());
K<s> := NumberField(x^2 - 2);

_<y> := PolynomialRing(K);
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


//_<z> := PolynomialRing(L);
//L2<ii,s7> := ext<L|z^2 + 1,z^2 + 7>;


// scott-Grassl 7a:
CC := ((-s+2)*r1+(3*s-2))*ii-s*r1-s;
AA := ((s+2)*r1+(s+4))*ii+(s-2)*r1-s;
BB := 4;

z8 := (1+ii)/s;

A := AA/CC;
B := BB/CC;

// scott-Grassl 7b"
CC2 := s*rr-s-2;
AA2 := 2;
BB2 := -rr+s-1;

A2 := AA2/CC2;
B2 := BB2/CC2;

psia := [1, A, A, B, A, B, B];
psib := [1,A2,A2,B2,A2,B2,B2];

u := (-r - 1 - s)/2;

C3 := 1;
A3 := u;
B3 := u^-1;

psib3 := [1,A3,A3,B3,A3,B3,B3];

renorm := func<psi | [psi[k]/psi[1] : k in [1..d]]>;

R2<a,b> := PolynomialRing(Rationals(),2);
phi2 := hom<R -> R2 | [x*y : x,y in [1,a,a,b,a,b,b]]>;
R2eqns := {@ phi2(f) : f in eqns @};

CC4 := 1/s7 + (1/s7 + s7)*(1+s)/2;
AA4 := (1+s)/(2*s7);
BB4 := -r/2;

A4 := AA4/CC4;
B4 := BB4/CC4;
	
check := func<v|[&+[v[1 + i]*v[1 + ((i + j) mod d)] : i in [0..d-1]] :  j in [0..d-1]]>;

OF := MaximalOrder(F);