

names := ["Y" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];
//Znames := ["Z" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];

Q := QuadraticField(2);
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

fp := func<i,j|f(i,j) + f(j,i)>;
fm := func<i,j|f(i,j) - f(j,i)>;

//S := Scheme(PZW, [f(i,j) : i in [0..d-1], j in [0..d-1] ]);
eqns := {@ g : g in {f(i,j): i in [0..d-1], j in [0..d-1] | [i,j] ne [0,0]} @};
eqns := [f : f in eqns];

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

RH := ideal<R|eqns>;
Idet := ideal<R|eqns2>;
Ra := ideal<R|a0-1>;
RY := ideal<R|Y(0,0)-1>;

Ia := RH + Idet + Ra;
IY := RH + Idet + RY;
aa := &+[Y(i+1,i) : i in [0..d-1]]; // Superconformal maps / d-symmetric spaces???

function  tt(a)
    t := ZeroMatrix(Q,d^2,d^2);
    for i in Integers(d) do
	for j in Integers(d) do
	    t[1 + d*(Integers() ! i) + (Integers() ! j), 1 + d*(Integers() ! (a*i)) + Integers() ! (a*j)] := 1;
	end for;
     end for;
    return t;
end function;

eqnsc := {@ Y(i,j) - Y(j,i) : i,j in [0..d] @};
Ic := ideal<R|eqnsc>;

eqnscc := {@ Y(i,j) - Y(-j,-i) : i,j in [0..d] @};
Icc := ideal<R|eqnscc>;

function IT(a)
    eqnsa := {@ Y(i,j) - Y(i,j)^tt(a) : i,j in Integers(d) @};
    return ideal<R|eqnsa>;
end function;


tc := ZeroMatrix(Q,d^2,d^2);
tw := ZeroMatrix(Q,d^2,d^2);
for i in Integers(d) do
    for j in Integers(d) do
	tc[1 + d*(Integers() ! i) + (Integers() ! j), 1 + d*(Integers() ! j) + Integers() ! i] := 1;
	tw[1 + d*(Integers() ! i) + (Integers() ! j), 1 + d*(Integers() ! j) + Integers() ! i] := 1;
    end for;
end for;

Tc := MatrixGroup<d^2,Q|tc>;

function TT(a)
    return MatrixGroup<d^2,Q|tt(a)>;
end function;
function TTc(a)
    return MatrixGroup<d^2,Q|tt(a),tc>;
end function;

//eqns3 := {@ Y(i,j) - Y(2*i,2*j) : i,j in Integers(d) @};

//I3 := ideal<R|eqns,eqns2,eqns3,a0-1>;
//I3 := ideal<R|eqns3>;


//T3 := MatrixGroup<d^2,Rationals()|t3>;

f0 := (R/Idet) ! (4*ff(0,0) - a0^2);
f1 := (R/Idet) ! (8*(ff(0,1) + ff(0,1) + ff(0,3)) - 3*a0^2);
f2 := (R/Idet) ! (ff(1,1) + ff(2,2) + ff(3,3));
f3 := (R/Idet) ! (ff(1,2) + ff(1,4) + ff(2,3));
f4 := (R/Idet) ! (ff(1,3) + ff(1,5) + ff(2,4));
f5 := (R/Idet) ! (ff(3,4) + ff(4,5) + ff(5,6));

fs := [f0,f1,f2,f3,f4,f5];

