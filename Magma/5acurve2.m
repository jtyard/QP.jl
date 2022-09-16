//Q3 := QuadraticField(3);
//R<x> := PolynomialRing(Q3);
//Kinf := AbsoluteField(ext<Q3 | x^2+1>);
//I := Kinf.1;
//s3 :=  Q3.1;
//w3 := I*s3;


d:=5;
K := QuadraticField((d-3)*(d+1));
OK := RingOfIntegers(K);

"Creating ray class group";
rayclassgroup, representing_ideal := RayClassGroup(d*OK,[1,2]);

"Building abelian extension";
FA := AbelianExtension(representing_ideal);

"Constructing number field";
F:= NumberField(FA);
//F := AbsoluteField(F);

//"Computing dth root of unity";
//zd := -Coefficient(Factorization(PolynomialRing(F) ! CyclotomicPolynomial(d))[1][1],0);

Kinf := CyclotomicField(12);
z12 := Kinf.1;                                                              
I := z12^3;
z3 := z12^4;
w3 := 1+2*z3;
s3 := I*w3;

OKinf := MaximalOrder(Kinf);
C12:= Kinf;
Z12:= OKinf;
C4 := sub<Kinf|I>;
Z4 := MaximalOrder(C4);
C3 := sub<Kinf|z3>;
Z3 := MaximalOrder(C3);
K := sub<Kinf|s3>;
OK := MaximalOrder(K);


//"Defining constants";
//s3 := Sqrt(F! 3);
//I := Sqrt(F!-1);
//w3 := I*s3;     

//"Constructing Kinf";
//Kinf := AbsoluteField(sub<F|w3,s3>);

// GKinf := GaloisGroup(Kinf);


g2 := (1125/2)*(3+4*I)*(1-w3);
g3 := 4125*s3*(11-2*I);

a := Kinf ! -g2/4;
b := Kinf ! -g3/4;

"Constructing elliptic curve";
E := EllipticCurve([a,b]);

//"Computing minimal model";
//Emin := MinimalModel(E);
jInvariant(E);

D := Discriminant(E);
C := Conductor(E);
Dfac := Factorization(ideal<OKinf|D>);
Cfac :=	Factorization(ideal<OKinf|C>);
//Factorization(ideal<OKinf|Discriminant(Emin)>);
//Factorization(ideal<OKinf|Conductor(Emin)>);


//f := DivisionPolynomial(E,5);
//fmin := DivisionPolynomial(Emin,5);

//"Factoring division polynomial of E";
//#Factorization(f);
//"Factoring division polynomial of Emin";
//#Factorization(fmin);

// Division polynomial factors completely over F(10^{1/3}) but not F.
//Fabs := AbsoluteField(F);
//"Factoring division polynomial of E over Kinf";
//#Factorization(PolynomialRing(Fabs) ! f);
//"Factoring division polynomial of Emin over Kinf";
//#Factorization(PolynomialRing(Fabs) ! fmin);

//R<x> := PolynomialRing(F);
//F10 := AbsoluteField(ext<F|x^3 - 10>);
//"Factoring division polynomial of E over Fi";
//#Factorization(PolynomialRing(F10) ! f);
//"Factoring division polynomial of Emin over Fi";
//#Factorization(PolynomialRing(F10) ! fmin);





P2 := Dfac[1,1];
P3 := Dfac[2,1];
P5 := Dfac[3,1];
PP5 := Dfac[4,1];
ideal<OKinf | D> eq P2^4*P3^12*P5^12*PP5^6;

//Discriminant(E) eq g2^3 - 27*g3^2;

//Emin := MinimalModel(E);

