
Q := Rationals();
Zd, ZZtoZd := Integers(d);

znames := ["z_" cat IntegerToString(i) : i in [0..d-1]];
wnames := ["w_" cat IntegerToString(i) : i in [0..d-1]];
Ynames := ["Y" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];

PZ<[z]> := ProjectiveSpace(Q,d-1);
AssignNames(~PZ,znames);

PW<[w]> := ProjectiveSpace(Q,d-1);
AssignNames(~PW,wnames);

PZW,phi := SegreProduct([PZ,PW]);
AssignNames(~PZW,Ynames);
R := CoordinateRing(PZW);


YY := [PZW.i : i in [1..d^2]];
Y := func<i,j| YY[d*(Integers() ! (Zd ! i)) + (Integers() ! (Zd ! j)) + 1]>;

//PZW := ProductProjectiveSpace(Q,[d-1,d-1]);
//AssignNames(~PZW, names);
zz := [PZ.i : i in [1..d]];
z := func<i| zz[(Integers() ! (Zd ! i)) + 1]>;

ww := [PW.i : i in [1..d]];
w := func<i| ww[(Integers() ! (Zd ! i)) + 1]>;

norm2 := &+[Y(i,i) : i in [0..d-1]];

function delta(i)
    if (Integers(d) ! i) eq 0 then
	return 1;
    else
	return 0;
    end if;
end function;

ff := func<i,j|&+[Y(k,k+i)*Y(k+i+j,k+j) : k in [0..d-1]]  >;

f := func<i,j| (d+1)*ff(i,j) - (delta(i) + delta(j))*norm2^2  >;

//S := Scheme(PZW, [f(i,j) : i in [0..d-1], j in [0..d-1] ]);

//ffeqns := {@ ff(i,j) : i in [0..d-1], j in [0..d-1] | i ge j @};
eqns := {@ f(i,j) : j in [-i..i], i in [0..(d+1)/2]  @};
S := Scheme(PZW, [f : f in eqns]);

RH := ideal<R|eqns>;
Ra := ideal<R|norm2 - 1>;

//time SS := Scheme(PZW, Radical(DefiningIdeal(S)));

function Del(f)
    return &+[Derivative(Derivative(f,z(k)),w(k)) : k in [0..d-1]];
end function;

//monomials := MonomialsOfDegree(R,4);

//function fvec(i,j)
//    fij := f(i,j);
//    return [MonomialCoefficient(fij,m) : m in monomials];
//end function;

//function ffvec(i,j)
//    ffij := ff(i,j);
//    return [MonomialCoefficient(ffij,m) : m in monomials];
//end function;

//N := #monomials;

//QN  := RModule(Integers(),N);

//GL2 := GL(2,Integers(d));
//SL2 := SL(2,Integers(d));


//Qd := CyclotomicField(d);
//function beta(i,j)
//    return &+[Qd.1^(j*(a-b))*z(a)*w(a+i)*z(b)*w(b+i) : a,b in [0..d-1]];
//end function;


function  tt(a)
    t := ZeroMatrix(Rationals(),2*d,2*d);
    for i in [0,1] do
	    for j in Integers(d) do
	        t[1 + d*i + (Integers() ! j), 1 + d*i + (Integers() ! (a*j))] := 1;
	    end for;    
    end for;
    return t;
end function;

eqnsc := {@ Y(i,j) - Y(j,i) : j in [0..i-1], i in [0..d] @};
Ic := ideal<R|eqnsc>;

eqnscc := {@ Y(i,j) - Y(-j,-i) : i,j in [0..d] @};
Icc := ideal<R|eqnscc>;

function IT(a)
    eqnsa := {@ Y(i,j)-Y(a*i,a*j) : i,j in [1..d-1]  @};
    return ideal<R|eqnsa>;
end function;

I := RH + IT(2) + Ic;
Ia := I + Ra;



tc := ZeroMatrix(Rationals(),2*d,2*d);
for i in [0..d-1] do
	tc[1 + d + i, 1 + i] := 1;
    tc[1 + i, 1 + d + i] := 1;
end for;

Tc := MatrixGroup<2*d,Rationals()|tc>;

function TT(a)
    return MatrixGroup<2*d,Rationals()|tt(a)>;
end function;
function TTc(a)
    return MatrixGroup<2*d,Rationals()|tt(a),tc>;
end function;