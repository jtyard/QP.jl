d:=5;
K := QuadraticField((d-3)*(d+1));
OK := RingOfIntegers(K);
"Creating ray class group";
rayclassgroup, representing_ideal := RayClassGroup(d*OK,[1,2]);
"Building abelian extension";
FA := AbelianExtension(representing_ideal);
"Constructing number field";
F:= NumberField(FA);

gal := func<s | ArtinMap(FA)(representing_ideal(s))>;


// Find (or define) a primitive dth roots of unity
// zd := F ! RootOfUnity(d,OptimizedRepresentation(AbsoluteField(F)));
"Finding dth root of unity";
zd := Coefficient(Factorization(PolynomialRing(F) ! CyclotomicPolynomial(d))[1][1],1);

"Computing constants";
// Define some constants
w3:=Sqrt(F!3);
w5:=Sqrt(F!5);
w15:= w3*w5;
I:=Sqrt(F!-1);

r1:=Sqrt(1/2*w5+5/2);
r2:=Sqrt((-1/8*w3+1/8*w5+3/8)*r1+1/8*w3*w5+5/8*w3+1/16*w5-5/16);

phi:=[
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


phi := Matrix(d,1,phi);

// Define some Galois group elements (at least complex conjugation)
c := gal(rayclassgroup.2);
sig := gal(rayclassgroup.1);
cc := sig^4*c;

// squared modulus
abs2:=func<x|x*c(x)>;

// inner product of two vectors
ip:=func<a,b|&+[A[i]*c(B[i]):i in [1..#Eltseq(a)]] where A:=Eltseq(a) where B:=Eltseq(b)>;

// Define the indexing group and module
sl2 := SL(2,Integers(d));
M := GModule(sl2);

// Define generalized Pauli matrices
X:=MatrixRing(F,d)!0;
X[1,d]:=1;
for i:=1 to d-1 do
    X[i+1,i]:=1;
end for;
Z := DiagonalMatrix(F,[zd^i:i in [0..d-1]]);

XZ := func<m| X^(Integers() ! m[1])*Z^(Integers() ! m[2]) >;

z2d := zd^(Integers() ! ((Integers(d) ! 2)^(-1)));
Delta := func<m| (-z2d)^(Integers() ! (m[1]*m[2]))*XZ(m) >;


phinorm := ip(phi,phi);
	
		   

f := map<M -> F | m :-> ip(phi,XZ(m)*phi)/phinorm >;
f2 := map<M -> F | m :-> ip(phi,Delta(m)*phi)/phinorm >;
act := func<g | map<M -> M | m :-> m * g> >;
		   


      

