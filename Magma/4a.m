d:=4;
dd:=8;

K := QuadraticField((d-3)*(d+1));
OK := RingOfIntegers(K);

"Creating ray class group";
rayclassgroup, representing_ideal := RayClassGroup(dd*OK,[1,2]);

"Building abelian extension";
FA := AbelianExtension(representing_ideal);

"Constructing number field";
F:= NumberField(FA);

gal := func<s | ArtinMap(FA)(representing_ideal(s))>;

Gal := [gal(s) : s in rayclassgroup];
c := Gal[3];

// Find (or define) a primitive dth roots of unity
"Computing ddth root of unity";
zdd := F ! RootOfUnity(2*d,OptimizedRepresentation(AbsoluteField(F)));

zd := zdd^2;

w2:= -Sqrt(F!2);
w5:= -Sqrt(F!5);
w10:=w2*w5;
r1:=Sqrt(w5+1);
I:=zd;

phiv:=[
8,
((w10+w2-2*w5-2)*r1-4)*I+(w10+w2)*r1+4*w2-4,
(8*w2-8)*I,
((w10+w2-2*w5-2)*r1+4)*I+(w10+w2)*r1-4*w2+4
];

phi := Matrix(d,1,phiv);
phiconj := Matrix(1,d,[c(a) : a in phiv]);

function apply_gal(M,s)
    for i in [1..Nrows(M)] do
	for j in [1..Ncols(M)] do
	    M[i,j] := s(M[i,j]);
	end for;
    end for;
    return M;
end function;


ctranspose := func<M | apply_gal(Transpose(M),c) >;

P0 := phi*phiconj; P0 := P0/Trace(P0);



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

z2d := -zdd;

Delta := func<j| (-z2d)^(Integers() ! (j[1]*j[2]))*XZ(j) >;

P := func<j|Delta(j)*P0*Delta([-j[1],-j[2]])>;



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

abs2(f2(M![1,1]));
