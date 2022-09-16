d:=3;



// Find (or define) a primitive dth roots of unity
// zd := F ! RootOfUnity(d,OptimizedRepresentation(AbsoluteField(F)));
"Computing dth root of unity";
//F := CyclotomicField(9);
//zdd := F.1;
//zd := zdd^3;

F := CyclotomicField(60);
zm := F.1;
c := Automorphisms(F)[4];
// c := Automorphisms(F)[2];                                                   // c := Automorphisms(F)[6];
zd := zm^20;
z4 := zm^15;
z3 := zd;
z5 := zm^12;
z12 := zm^5;

//phi:=[F ! 0,F ! 1,F ! (-1)];
al := F! z12^-1;
phi := [F!0,al^-1,-al];
phi := Matrix(d,1,phi);


"Defining other functions and operators";

// Define some Galois group elements (at least complex conjugation)

// squared modulus
abs2:=func<x|x*c(x)>;

// inner product of two vectors
ip:=func<a,b|&+[A[i]*c(B[i]):i in [1..#Eltseq(a)]] where A:=Eltseq(a) where B:=Eltseq(b)>;

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

phinorm := ip(phi,phi);
	
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

// Coset representatives for the map s -> gl2/S
// Don't need anything else because cc acts trivially
// Therefore c acts as sig^4.

phis := {};
for m in M do
    phis := phis join {Delta(m)*phi};
end for;

V := VectorSpace(F,d);

ranks := [Rank(sub<V | [V ! Eltseq(phi) : phi in subs]  >) : subs in Subsets(phis,3)];

S3 := [s : s in Subsets({m : m in M},3) | Rank(sub<V | [V ! Eltseq(Delta(m)*phi) : m in s]>) eq 2];
#S3;

function count(A)
    S := {x : x in A};
    c := {};
    for x in S do
	c := c join {[x,#[y : y in A | x eq y]]};
    end for;
    return c;
end function;

count(ranks);
	
		   




      

