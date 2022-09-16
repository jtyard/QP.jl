d:=3;

F := CyclotomicField(60);
//zm := F.1;
// c := Automorphisms(F)[4];
// c := Automorphisms(F)[2];
// c := Automorphisms(F)[6];
c := ComplexConjugate;

zm := F.1;
zd := zm^20;
z4 := zm^15;
z3 := zd;
z5 := zm^12;
z12 := zm^5;
z6 := z12^2;

//philist := Eltseq(Vector([F!0, 1,-al])*Matrix([[z3^(i*j) : j in [0..2] ] : i in [0..2]]));

philist := [F!0, 1,al];


"Defining other functions and operators";

// Define some Galois group elements (at least complex conjugation)

// squared modulus
abs2:=func<x|x*c(x)>;

// inner product of two vectors
ip:=func<a,b|&+[c(A[i])*B[i]:i in [1..#Eltseq(a)]] where A:=Eltseq(a) where B:=Eltseq(b)>;

// Define the indexing group and module
sl2 := SL(2,Integers(d));
gl2 := GL(2,Integers(d));
M := GModule(gl2);

// Define generalized Pauli matrices
X:=MatrixRing(F,d)!0;
X[1,d]:=1;
for i:=1 to d-1 do
    X[i+1,i]:=1;
end for;
Z := DiagonalMatrix(F,[zd^i:i in [0..d-1]]);

XZ := func<j| X^(Integers() ! j[1])*Z^(Integers() ! j[2]) >;

z2d := -zd^(Integers() ! ((Integers(d) ! 2)^(-1)));

Delta := func<j| (-z2d)^(Integers() ! (j[1]*j[2]))*XZ(j) >;

phinorm := ip(philist,philist);

phi := Matrix(d,1,philist);

"Defining other functions and operators";

function apply_gal(M,s)
    for i in [1..Nrows(M)] do
	for j in [1..Ncols(M)] do
	    M[i,j] := s(M[i,j]);
	end for;
    end for;
    return M;
end function;
	    
ctranspose := func<M | apply_gal(Transpose(M),c) >;
		   
P0 := phi*ctranspose(phi)/phinorm;

phiconj := Vector([c(x) : x in philist]);

sic := [[Delta([a,b])*phi : b in [0..d-1]] : a in [0..d-1]];

overlaps := [[phiconj*sic[a+1][b+1]/phinorm : b in [0..d-1]] : a in [0..d-1]];

P := [[Delta([a,b])*P0*Delta([-a,-b]) : b in [0..d-1]] : a in [0..d-1]];

function Tr3(i,j,k)
    return Trace(P[i[1]+1,i[2]+1]*P[j[1]+1,j[2]+1]*P[k[1]+1,k[2]+1]);
end function;


"Computing triple phases";
phases := [[[[ Tr3([0,0], [j1,j2], [k1,k2]) : k2 in [0..d-1]] : k1 in [0..d-1]] : j2 in [0..d-1]] : j1 in [0..d-1]];

phaseset := {@ phases[j1,j2,k1,k2] : j1,j2,k1,k2 in [1..d] | [j1,j2] ne [1,1] and [k1,k2] ne [1,1] and [j1,j2] ne [k1,k2] @};

//pstemp := phaseset;
//sigsc := [c^j*sig^k : k in [0..7], j in [0,1]];
//orbits := [];
//while pstemp ne {} do
//    a := pstemp[1];
//    orbit := {@ sigsc[k](a) : k in [1..16] @};
//    Append(~orbits,orbit);
//    pstemp := pstemp diff orbit;
//end while;
//
//sigs := [sig^k : k in [0..7]];
//realset := {@ (a+c(a))/2 : a in phaseset @};
//rstemp := realset;
//rorbits := [];
//while rstemp ne {} do
//    r:= rstemp[1];
//    orbit := {@ sigs[k](r) : k in [1..8] @};
//    Append(~rorbits,orbit);
//    rstemp := rstemp diff orbit;
//end while;

function phasetopairs(a)
    return [[M ! [j1,j2],M! [k1,k2]] : j1, j2, k1, k2 in [0..d-1] | phases[j1+1][j2+1][k1+1][k2+1] eq a];
end function;

if Min([Trace(P0*P[j,k]) eq 1/4 : j,k in [1..d] | (j ne 1) or (k ne 1)]) eq true then
    "Is a SIC-POVM";
else
    "NOT a SIC-POVM";
end if;

"Number of triple products for each root of unity";

roots := [zm^k : k in [0..59]];

phasepowers := [Index(roots,8*a)-1 : a in phaseset];
Sort(~phasepowers);
trips := [[k/60,#phasetopairs(zm^k/8)] : k in phasepowers];
trips;
&+[a[2] : a in trips];


// Compute ranks of spans of triples

phis := {};
for m in M do
    phis := phis join {Delta(m)*phi};
end for;

V := VectorSpace(F,d);

ranks := [Rank(sub<V | [V ! Eltseq(phi) : phi in subs]  >) : subs in Subsets(phis,3)];

function count(A)
    S := {x : x in A};
    c := {};
    for x in S do
	c := c join {[x,#[y : y in A | x eq y]]};
    end for;
    return c;
end function;

"Number of triangles of a given rank";
count(ranks);
	
		   
f := map<M -> F | m :-> ip(phi,XZ(m)*phi)/phinorm >;
f2 := map<M -> F | m :-> ip(phi,Delta(m)*phi)/phinorm >;
act := func<g | map<M -> M | m :-> m * g> >;

function equal(f,g)
    for m in M do
	if f(m) ne g(m) then
	    return false;
	end if;
    end for;
    return true;
end function;

// Stabilizer group
stab := [g : g in gl2 | equal(f2,act(g)*f2)];
"|gl2 stabilizer|=", #stab;

slstab := [g : g in stab | Determinant(g) eq 1];
elstab := slstab cat [g : g in stab | Determinant(g) eq -1];

// Coset representatives for the map s -> gl2/S
// Don't need anything else because cc acts trivially
// Therefore c acts as sig^4.

// Not the usual action.  Applies m first then applies g.
function gramact(gm,gram)
    out := ZeroMatrix(F,d,d);
    g := gm[1];
    m := gm[2];
    for k in [0..2] do
	for j in [0..2] do
	    mm := (M! [j,k])*Transpose(g^-1);
	    a := Integers() ! mm[1];
	    b := Integers() ! mm[2];
	    out[j+1,k+1] := z3^((Integers() ! m[2])*j - (Integers() ! m[1])*k)*gram[a + 1, b + 1];
	end for;
    end for;
    return out;
end function;

mphase := func<m,m0|z3^(Integers() ! (m0[2]*m[1] - m0[1]*m[2]))>;

gram := Matrix([[f(M ! [j,k]) : k in [0..d-1]] : j in [0..d-1]]);

Matrix([[#[g : g in gl2 | gramact(<g,M![j,k]>,gram) eq gramact(<gl2!1,M![j,k]>,gram)] : k in [0,1,2]] : j in [0,1,2]]);





      





      


    
		   




      

