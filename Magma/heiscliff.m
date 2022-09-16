// d:=3;
Qd := CyclotomicField(d);
ZQd := MaximalOrder(Qd);

zd := Qd.1;

V := VectorSpace(Qd,d);
A := Algebra(MatrixAlgebra(Qd,d));

// Define the indexing group and module
sl2 := SL(2,Integers(d));
gl2 := GL(2,Integers(d));
M := GModule(gl2);

// Define generalized Pauli matrices
X:=MatrixRing(Qd,d)!0;
X[1,d]:=1;
for i:=1 to d-1 do
    X[i+1,i]:=1;
end for;
Z := DiagonalMatrix(Qd,[zd^i:i in [0..d-1]]);

XZ := func<j| X^(Integers() ! j[1])*Z^(Integers() ! j[2]) >;

z2d := -zd^(Integers() ! ((Integers(d) ! 2)^(-1)));

Delta := func<j| (-z2d)^(Integers() ! (j[1]*j[2]))*XZ(j) >;

hc := func<j,k|(-z2d)^(Integers() ! (j[2]*k[1] - j[1]*k[2])) >;

relations := [<d*(Integers()!j[1]) + (Integers()!j[2])+1,d*(Integers()!k[1]) + (Integers()!k[2])+1,d*(Integers()!(j+k)[1]) + (Integers()!(j+k)[2]) +1, hc(j,k)> : k in M, j in M];

//A := AssociativeAlgebra<Qd,d^2 |relations>;
Amat := MatrixAlgebra(Qd,d);
A, f := Algebra(Amat);
//MA := MaximalOrder(A);
//MB := MaximalOrder(B);

O := Order([f(Amat ! Delta(j)) : j in M]);
MO := MaximalOrder(O);
DO := Discriminant(O);
DMO := Discriminant(MO);
DO;
DMO;
d^(d^2), d^(2*d^2);

// So MO <-> Hurwitz quaternions, with O <-> Lipschitz...

H:= MatrixGroup<d,Qd|X,Z>;

B := ZBasis(O);
Gram := Matrix(RealField(), [[AbsoluteTrace(Trace((f^-1)(b1)*(f^-1)(b2))) : b2 in B] : b1 in B]);

BM := ZBasis(MO);
GramM := Matrix(RealField(), [[AbsoluteTrace(Trace((f^-1)(b1)*(f^-1)(b2))) : b2 in BM] : b1 in BM]);


NumericalEigenvalues(Matrix(RealField(),Gram));
NumericalEigenvalues(Matrix(RealField(), GramM));
