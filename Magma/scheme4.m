
znames := ["z" cat IntegerToString(i) : i in [0..d-2]];
wnames := ["w" cat IntegerToString(i) : i in [0..d-2]];

names := znames cat wnames;
Q := Rationals();
// PZ<[z]> := ProjectiveSpace(Q,d-1);
// PW<[w]> := ProjectiveSpace(Q,d-1);
// PZW := DirectProduct(PZ,PW);

AZW := AffineSpace(Q,2*d-2);
AssignNames(~AZW, names);
zz := [AZW.i : i in [1..d-1]] cat [1];
ww := [AZW.i : i in [d..2*d-2]] cat [1];

Zd, ZZtoZd := Integers(d);

z := func<i| zz[(Integers() ! (Zd ! i)) + 1]>;
w := func<i| ww[(Integers() ! (Zd ! i)) + 1]>;

norm2 := &+[z(i)*w(i) : i in [0..d-1]];

function delta(i)
    if (Integers(d) ! i) eq 0 then
	return 1;
    else
	return 0;
    end if;
end function;

ff := func<i,j|&+[z(k)*w(k+i)*w(k+j)*z(k+i+j) : k in [0..d-1]]  >;

f := func<i,j| (d+1)*ff(i,j) - (delta(i) + delta(j))*norm2^2  >;

//S := Scheme(PZW, [f(i,j) : i in [0..d-1], j in [0..d-1] ]);

ffeqns := {ff(i,j) : i in [0..d-1], j in [0..d-1] | i ge j};
eqns := {@ f(i,j) : i in [0..d-1], j in [0..d-1] | i ge j @};
S := Scheme(AZW, [f : f in eqns]);

GL2 := GL(2,Integers(d));
SL2 := SL(2,Integers(d));

