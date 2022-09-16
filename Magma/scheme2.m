

names := ["Y" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];
//Znames := ["Z" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];

Q := Rationals();
// PZ<[z]> := ProjectiveSpace(Q,d-1);
// PW<[w]> := ProjectiveSpace(Q,d-1);
// PZW := DirectProduct(PZ,PW);

PY := ProjectiveSpace(Q,d^2-1);
AssignNames(~PY, names);
YY := [PY.i : i in [1..d^2]];
//ww := [PZW.i : i in [d+1..2*d]];


Zd, ZZtoZd := Integers(d);

Y := func<i,j| YY[d*(Integers() ! (Zd ! i)) + (Integers() ! (Zd ! j)) + 1]>;
//w := func<i| ww[(Integers() ! (Zd ! i)) + 1]>;

norm2 := &+[Y(i,i) : i in [0..d-1]];

function delta(i)
    if (Integers(d) ! i) eq 0 then
	return 1;
    else
	return 0;
    end if;
end function;

ff := func<i,j|&+[Y(k,k+i)*Y(k+i+j,k+j) : k in [0..d-1]] >;
f:=func<i,j|(d+1)*ff(i,j) - (delta(i) + delta(j))*norm2^2  >;

//S := Scheme(PZW, [f(i,j) : i in [0..d-1], j in [0..d-1] ]);
eqns := [g : g in {f(i,j): i in [0..d-1], j in [0..d-1] }];

eqns2 := [Y(i,k)*Y(j,l) - Y(i,l)*Y(j,k) : i,j,k,l in [0..d-1] | i lt j and k lt l];

//eqns3 := [norm2^2 - &+[Y(i,j)*Y(j,i) : i,j in [0..d-1]]];

S := Scheme(PY, eqns cat eqns2);

hc := func<p | [ComplexConjugate(p[j*d + i + 1]) : i, j in [0..d-1]]>;

function FtoM(p)
    m := Matrix([[p[d*i + j + 1] : j in [0..d-1]] : i in [0..d-1]]);
    return m/Trace(m);
end function;

//time SS := Scheme(PZW, Radical(DefiningIdeal(S)));

//function Del(f)
//    return &+[Derivative(Derivative(f,z(k)),w(k)) : k in [0..d-1]];
//end function;

//carpow := CartesianPower([0,1,2],d);
//indices := [[A[k] : k in [1..d]]  : A in carpow | &+[A[k] : k in [1..d]] eq 2];
//monomials := [&*[z(k)^A[k+1] : k in [0..d-1]] * &*[w(k)^B[k+1] : k in [0..d-1]]  : A,B in indices];

//function fvec(i,j)
//    ff := f(i,j);
//    return [MonomialCoefficient(ff,m) : m in monomials];
//end function;

//N := #monomials;

//QN := RModule(Rationals(),N);
//M:= sub<QN|[QN ! fvec(i,j) : i,j in [0..d-1] | i le j]>;

if IsEven(d) then
    harminv := ((d+3)*(d-1) + 3)/4;
else
    harminv := (d+3)*(d-1)/4;
end if;

//eqns2 := [ &+[(QN ! b)[i]*monomials[i] : i in [1..N]] : b in Basis(M)];

//S2 := Scheme(PZW,eqns2);

print harminv;

sl2 := SL(2,Integers(d));


a0 := &+[Y(i,i) : i in [0..d-1]];

R := CoordinateRing(PY);
AY := AffineSpace(R);

I := ideal<R|eqns,eqns2>;
Ia := ideal<R|I,a0-1>;

a := &+[Y(i+1,i) : i in [0..d-1]]; // Superconformal maps / d-symmetric spaces???

eqns3 := {@ Y(i,j) - Y(2*i,4*j) : i,j in Integers(7) @};

//I3 := ideal<R|eqns,eqns2,eqns3,a0-1>;
I3 := ideal<R|eqns3>;

eqnsgood1 := {@ Y(i,j) - Y(5*j,3*i) : i,j in [0..6] @};
eqnsgood2 := {@ Y(i,j) - Y(2*j,4*j) : i,j in [0..6] @};

eqnsc    := {@ Y(i,j) - Y(j,i) : i,j in [0..6] @};
eqnscc   := {@ Y(i,j) - Y(-j,-i) : i,j in [0..6]  @};
eqns6    := {@ Y(i,j) - Y(3*i,5*j) : i,j in [0..6] @};
eqnsccc  := {@ Y(i,j) - Y(-i,-j) : i,j in [0..6] @};
//I6 := ideal<R|eqns,eqns2,eqns6,a0-1>;
//I3c := ideal<R|eqns,eqns2,eqns3,eqnsc,a0-1>;
//Ic := ideal<R|eqns,eqns2,eqnsc,a0-1>;
//Icc := ideal<R|eqns,eqns2,eqnscc,a0-1>;

Igood1 := ideal<R|eqnsgood1>;
Igood2 := ideal<R|eqnsgood2>;
Ic := ideal<R|eqnsc>;
Icc := ideal<R|eqnscc>;
I6 := ideal<R|eqns6>;
Iccc := ideal<R|eqnsccc>;
I12 := I6 + Ic;
I12 eq I6 + Icc;

// IsPrime(Ia + I3 + Ic);
// IsPrime(Ia + I3 + Icc);
eqnsj3 := {@ Y(i,j) - Y(2*i,2*j) : i,j in [0..6] @};
J3 := ideal<R|eqnsj3>;
J3a := ideal<R|eqnsj3,a0-1>;




t3 := Matrix(Rationals(),d^2,[0 : i in [1..d^4]]);
for i in Integers(d) do
    for j in Integers(d) do
	t3[1 + d*(Integers() ! i) + Integers() ! j, 1 + d*(Integers() ! (2*i)) + Integers() ! (2*j)] := 1;
    end for;
end for;

T3 := MatrixGroup<d^2,Rationals()|t3>;

