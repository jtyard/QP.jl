// 3-23-20
// Build the Scott-Grassl 7a and 7b fiducials once and for all
// Confusion: For 7b to be real, r1 i needs to be real, which means that w2 = -Sqrt(2)
// So I should read r1 i as Sqrt(- (-2 s + 1)) = Sqrt(2s - 1)

d:=7;

_<xx> := PolynomialRing(Rationals());
K<s> := NumberField(xx^2 - 2);

Kx<x> := PolynomialRing(K);
F<r,ii> := ext<K|x^2 - (2*s-1), x^2 + 1>;
Fabs := AbsoluteField(F);

s := Fabs ! s;
r := Fabs ! r;
ii := Fabs ! r;
//s7 := 1 + 2*(zd + zd^2 + zd^4);

// Here, my r is SG's r1*I, so I also define
r1 := ii*r;

Gal := Automorphisms(Fabs);

function theaut(gr,gii)
    for g in Gal do
        if g(r) eq gr and g(ii) eq gii then
            return g;
        end if;
    end for;
end function;

c1 := theaut(r,-ii);
c2 := theaut(-r,-ii);
c := c1*c2;

//c1 := hom<F->F|r,-ii,zd^-1>;
//c2 := hom<F->F|-r,-ii,zd^-1>;

/*
// Define generalized Pauli matrices
Xd:=MatrixRing(F,d)!0;
Xd[1,d]:=1;
for i:=1 to d-1 do
    Xd[i+1,i]:=1;
end for;
Zd := DiagonalMatrix(F,[zd^i:i in [0..d-1]]);

XZ := func<j| Xd^(Integers() ! j[1])*Zd^(Integers() ! j[2]) >;

z2d := -zd^(Integers() ! ((Integers(d) ! 2)^(-1)));

Delta := func<j| (-z2d)^(Integers() ! (j[1]*j[2]))*XZ(j) >;

*/

//_<z> := PolynomialRing(L);
//L2<ii,s7> := ext<L|z^2 + 1,z^2 + 7>;


// Scott-Grassl 7a?

aa0 := ((-s+2)*r+(3*s-2))*ii-s*r-s;
aa1 := ((s+2)*r+(s+4))*ii+(s-2)*r-s;
aa3 := 4;

z8 := (1+ii)/s;

a1 := aa1/aa0;
a3 := aa3/aa0;

// Scott-Grassl 7b?
bb0 := s*r-s-2;
bb1 := 2;
bb3 := -r+s-1;

b1 := bb1/bb0;
b3 := bb3/bb0;

psia := [1, a1, a1, a3, a1, a3, a3];
psib := [1, b1, b1, b3, b1, b3, b3];

u := (r - 1 - s)/2;

C1 := u;
C3 := u^-1;

psib3 := [1,C1,C1,C3,C1,C3,C3];

renorm := func<psi | [psi[k]/psi[1] : k in [1..d]]>;

//R2<a,b> := PolynomialRing(Rationals(),2);
//phi2 := hom<R -> R2 | [x*y : x,y in [1,a,a,b,a,b,b]]>;
//R2eqns := {@ phi2(f) : f in eqns @};

//CC4 := 1/s7 + (1/s7 + s7)*(1+s)/2;
//AA4 := (1+s)/(2*s7);
//BB4 := -r/2;

//A4 := AA4/CC4;
//B4 := BB4/CC4;
	
/*
braket := func<w,v | &+[c1(w[i])*v[i] : i in [1..d]]>;
abs2 := func<a|a*c1(a)>;

function check(v)
    norm2 := braket(v,v);
    vvec := Vector(v);
    al := [[abs2(&+[v[1+a]*v[1+((a+j1) mod d)]*zd^(a*j2) : a in [0..d-1]]/norm2) :  j2 in [0..d-1]] : j1 in [0..d-1]];
    return Matrix(al);
end function;
*/

//OF := MaximalOrder(F);