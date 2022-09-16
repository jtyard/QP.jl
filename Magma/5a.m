d:=5;
K := NumberField(QuadraticField((d-3)*(d+1)));
OK := RingOfIntegers(K);

"Creating ray class group";
rayclassgroup, representing_ideal := RayClassGroup(d*OK,[1,2]);

"Building abelian extension";
FA := AbelianExtension(representing_ideal);

"Constructing number field";
F:= NumberField(FA);

V := VectorSpace(F,d);

art := func<s | ArtinMap(FA)(representing_ideal(s))>;

Gal := [art(s) : s in rayclassgroup];

// Define some Galois group elements (at least complex conjugation)           
c   := art(rayclassgroup.2);
sig := art(rayclassgroup.1);
cc  := sig^4*c;
Sig := [sig^i : i in [1..8]];

// squared modulus                                                           
abs2:=func<x|x*c(x)>;

// inner product of two vectors                                              
ip:=func<a,b|&+[A[i]*c(B[i]):i in [1..#Eltseq(a)]] where A:=Eltseq(a) where B:=Eltseq(b)>;


// Find (or define) a primitive dth roots of unity
// zd := F ! RootOfUnity(d,OptimizedRepresentation(AbsoluteField(F)));
"Computing dth root of unity";
zd := -Coefficient(Factorization(PolynomialRing(F) ! CyclotomicPolynomial(d))[1][1],0);

"Loading fiducial";
// Define some constants
w3:=Sqrt(F!3);
w5:=Sqrt(F!5);
w15:= w3*w5;
I:=Sqrt(F!-1);

r1:=Sqrt(1/2*w5+5/2);
r2:=Sqrt((-1/8*w3+1/8*w5+3/8)*r1+1/8*w3*w5+5/8*w3+1/16*w5-5/16);

philist:= [
    16*r1,
    (((8*w3-10*w5-2)*r1+20*w3-16*w5)*r2+((2*w5+2)*r1+(2*w5+10)))*I
        +((-6*w15-6*w3+8*w5+24)*r1+(-12*w15-20*w3+12*w5+40))*r2+(4*w5-8)*r1-4*w15-4*w5,
    (((16*w3-14*w5+6)*r1+(-2*w15+30*w3-24*w5))*r2+((w15-5*w3+4)*r1-3*w5+5))*I
        +((2*w15-2*w3+4*w5-12)*r1-2*w5-10)*r2+(w5+3)*r1+3*w15-5*w3-2*w5+10,
    (((-6*w15-6*w3+6*w5+14)*r1+(-8*w15-20*w3+16*w5+20))*r2+((w15+5*w3-2*w5+6)*r1+4*w15-2*w5))*I
        +((2*w15+2*w3-2*w5-18)*r1+20*w3-12*w5)*r2+(3*w5-1)*r1-2*w15,
    (((-8*w15+12*w3-12*w5+32)*r1+(-10*w15+10*w3-6*w5+30))*r2+((-2*w15+w5+3)*r1+w15+5*w3+3*w5+5))*I
        +((-4*w3-4)*r1+2*w15-10*w3+2*w5-10)*r2+(-w15+5*w3+2*w5)*r1+3*w15+5*w3+w5+5
];

phi := Matrix(d,1,philist);

"Defining other functions and operators";

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

phi := XZ([3,1])*phi;


phinorm := ip(phi,phi);

function apply_gal(M,s)
    for i in [1..Nrows(M)] do
	for j in [1..Ncols(M)] do
	    M[i,j] := s(M[i,j]);
	end for;
    end for;
    return M;
end function;
	    
ctranspose := func<M | apply_gal(Transpose(M),c) >;

		   

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
S := [ScalarMatrix(IntegerRing(5), 2, 1),
      Matrix(IntegerRing(5), 2, 2, [ 0, 2, 2, 4 ]),
      Matrix(IntegerRing(5), 2, 2, [ 4, 3, 3, 0 ])];

// Coset representatives for the map s -> gl2/S
// Don't need anything else because cc acts trivially
// Therefore c acts as sig^4.
Ssigs := [ScalarMatrix(IntegerRing(5), 2, 1),
	  Matrix(IntegerRing(5), 2, 2, [ 3, 3, 3, 4 ]),
	  ScalarMatrix(IntegerRing(5), 2, 2),
	  Matrix(IntegerRing(5), 2, 2, [ 2, 3, 3, 3 ]),
	  ScalarMatrix(IntegerRing(5), 2, 4),
	  Matrix(IntegerRing(5), 2, 2, [ 2, 2, 2, 1 ]),
	  ScalarMatrix(IntegerRing(5), 2, 3),
	  Matrix(IntegerRing(5), 2, 2, [ 3, 4, 4, 1 ])];

//P0 := Matrix([[philist[i]*c(philist[j]) : j in [1..d] ] : i in [1..d]]);
P0 := phi*ctranspose(phi)/phinorm;

phiconj := Vector([c(x) : x in philist]);

sic := [[Delta([a,b])*phi : b in [0..d-1]] : a in [0..d-1]];

overlaps := [[phiconj*sic[a+1][b+1]/phinorm : b in [0..d-1]] : a in [0..d-1]];

P := [[Delta([a,b])*P0*Delta([-a,-b]) : b in [0..d-1]] : a in [0..d-1]];

function Tr3(i,j,k)
    return Trace(P[i[1]+1,i[2]+1]*P[j[1]+1,j[2]+1]*P[k[1]+1,k[2]+1]);
end function;

al := zd^2*f2([1,0]);
