
names := ["Y" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];
//Znames := ["Z" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];

Q := Rationals();
// PZ<[z]> := ProjectiveSpace(Q,d-1);
// PW<[w]> := ProjectiveSpace(Q,d-1);
// PZW := DirectProduct(PZ,PW);

AY := AffineSpace(Q,d^2);
AssignNames(~AY, names);
YY := [AY.i : i in [1..d^2]];

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
f:=func<i,j|(d+1)*ff(i,j) - (delta(i) + delta(j))  >;

//S := Scheme(PZW, [f(i,j) : i in [0..d-1], j in [0..d-1] ]);
eqns := [g : g in {f(i,j): i in [0..d-1], j in [0..d-1] }];

minors := [Y(i,k)*Y(j,l) - Y(i,l)*Y(j,k) : i,j,k,l in [0..d-1] | i lt j and k lt l];

eqns3 := [g : g in {f(i,j) + f(j,i) : i,j in [0..d-1]}];

eqns4 := [g : g in {f(i,j) - f(j,i) : i,j in [0..d-1]}];
//eqns3 := [norm2^2 - &+[Y(i,j)*Y(j,i) : i,j in [0..d-1]]];

TrY := norm2;
TrY2 := &+[Y(i,j)*Y(j,i) : i,j in [0..d-1]];
TrY3 := &+[Y(i,j)*Y(j,k)*Y(k,i) : i,j,k in [0..d-1]];

// This one seems to be reduced
S := Scheme(AY,eqns cat [TrY2 - 1,TrY-1]);
S2 := Scheme(AY,eqns cat minors cat [TrY3-1, TrY2 - 1, TrY-1]);
S3 := Scheme(AY, eqns cat eqns3 cat [TrY^2 - 1, TrY2 - 1]);
S4 := Scheme(AY, eqns cat eqns4 cat [TrY^2 - 1, TrY2 - 1]);
S34 := Scheme(AY, eqns cat eqns3 cat eqns4 cat [TrY^2 - 1, TrY2 - 1]);

hc := func<p | [ComplexConjugate(p[j*d + i + 1]) : i, j in [0..d-1]]>;

function FtoM(p)
    m := Matrix([[p[d*i + j + 1] : j in [0..d-1]] : i in [0..d-1]]);
    return m/Trace(m);
end function;

