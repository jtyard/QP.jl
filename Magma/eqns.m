// Only define all the eqns in this file.

Q := Rationals();
Zd, ZZtoZd := Integers(d);

RR := PolynomialRing(Q,d^2);
Ynames := ["Y" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];
AssignNames(~RR, Ynames);
YY := [RR.i : i in [1..d^2]];
Y := func<i,j| YY[d*(Integers() ! (Zd ! i)) + (Integers() ! (Zd ! j)) + 1]>;

R := PolynomialRing(Q,2*d);
znames := ["z_" cat IntegerToString(i) : i in [0..d-1]];
wnames := ["w_" cat IntegerToString(i) : i in [0..d-1]];
zwnames := znames cat wnames;
AssignNames(~R,zwnames);
zz := [R.i : i in [1..d]];
ww := [R.i : i in [d+1..2*d]];
z := func<i| zz[(Integers() ! (Zd ! i)) + 1]>;
w := func<i| ww[(Integers() ! (Zd ! i)) + 1]>;

toR := hom<RR->R|[z(i)*w(j) : i,j in [0..d-1]]>;

TrX := &+[Y(i,i) : i in [0..d-1]];
TrX2 := &+[Y(i,j)*Y(j,i) : i,j in [0..d-1]];

A0 := TrX;
B0p := (TrX^2 + TrX2)/2;
B0m := (TrX^2 - TrX2)/2;

function delta(i)
    if (Integers(d) ! i) eq 0 then
	return 1;
    else
	return 0;
    end if;
end function;


ff := func<i,j|&+[Y(k,k+i)*Y(k+i+j,k+j) : k in [0..d-1]] >;
//f  := func<i,j|(d+1)*ff(i,j) - (delta(i) + delta(j))*TrX^2 >;
f  := func<i,j|(d^2-1)*ff(i,j) - (d*delta(i) - delta(j))*TrX^2 - (-delta(i) + d*delta(j))*TrX2 >;

ffp := func<i,j|(ff(i,j) + ff(j,i))/2>;
ffm := func<i,j|(ff(i,j) - ff(j,i))/2>;

fp := func<i,j|ffp(i,j) - (delta(i) + delta(j))*B0p/(d+1)>;
fm := func<i,j|ffm(i,j) - (delta(i) - delta(j))*B0m/(d-1)>;

eqns  := {@ f(i,j)  : j in [0..d-1], i in [0..d-1] @};
eqnsp := {@ fp(i,j) : j in [i..d-1], i in [0..d-1] @};
eqnsm := {@ fm(i,j) : j in [i..d-1], i in [0..d-1] @};

if IsEven(d) then
    else
    eqns2 := [f(i,j): j in [-i..i], i in [1..(d-1)/2]];
end if;

eqnsdet := [ Y(i,k)*Y(j,l) - Y(i,l)*Y(j,k) : i,j,k,l in [0..d-1] | i lt j and k lt l ];

RRH := ideal<RR|eqns>;
RRa := ideal<RR|A0-1>;
RRY := ideal<RR|Y(0,0)-1>;

Idet := ideal<RR|eqnsdet>;

Del := func<f|&+[Derivative(f,Y(i,i)) : i in [0..d-1]]>;


//S := Scheme(PY, eqns cat eqnsdet);

hc := func<p | [ComplexConjugate(p[j*d + i + 1]) : i, j in [0..d-1]]>;

function FtoM(p)
    m := Matrix([[p[d*i + j + 1] : j in [0..d-1]] : i in [0..d-1]]);
    return m/Trace(m);
end function;



sl2 := SL(2,Integers(d));

a0 := &+[Y(i,i) : i in [0..d-1]];

aa := &+[Y(i+1,i) : i in [0..d-1]]; // Superconformal maps / d-symmetric spaces???

function  tt(a)
    t := ZeroMatrix(Rationals(),d^2,d^2);
    for i in Integers(d) do
	for j in Integers(d) do
	    t[1 + d*(Integers() ! i) + (Integers() ! j), 1 + d*(Integers() ! (a*i)) + Integers() ! (a*j)] := 1;
	end for;
     end for;
    return t;
end function;

eqnsc := {@ Y(i,j) - Y(j,i) : i,j in [0..d] @};
Ic := ideal<RR|eqnsc>;

eqnscc := {@ Y(i,j) - Y(-j,-i) : i,j in [0..d] @};
Icc := ideal<RR|eqnscc>;

function IT(a)
    eqnsa := {@ Y(i,j) - Y(i,j)^tt(a) : i,j in Integers(d) @};
    return ideal<RR|eqnsa>;
end function;

tc := ZeroMatrix(Rationals(),d^2,d^2);
tw := ZeroMatrix(Rationals(),d^2,d^2);
for i in Integers(d) do
    for j in Integers(d) do
	tc[1 + d*(Integers() ! i) + (Integers() ! j), 1 + d*(Integers() ! j) + Integers() ! i] := 1;
	tw[1 + d*(Integers() ! -i) + (Integers() ! -j), 1 + d*(Integers() ! i) + Integers() ! j] := 1;
    end for;
end for;

Tc := MatrixGroup<d^2,Rationals()|tc>;

function TT(a)
    return MatrixGroup<d^2,Rationals()|tt(a)>;
end function;
function TTc(a)
    return MatrixGroup<d^2,Rationals()|tt(a),tc>;
end function;


if d eq 7 then
    T3c := MatrixGroup<d^2,Rationals()|tt(2),tc>;
    T6c := MatrixGroup<d^2,Rationals()|tt(3),tc>;
end if;

//eqns3 := {@ Y(i,j) - Y(2*i,2*j) : i,j in Integers(d) @};

//I3 := ideal<R|eqns,eqns2,eqns3,a0-1>;
//I3 := ideal<R|eqns3>;
//T3 := MatrixGroup<d^2,Rationals()|t3>;

monomialsY2 := MonomialsOfDegree(RR,2);
RRV := VectorSpace(Q,#monomialsY2);
RRv := func<f|RRV ! [MonomialCoefficient(f,m) : m in monomialsY2]>;
RRspan := func<fs|sub<RRV| [RRv(f) : f in fs]>>;
RRdim := func<fs|Dimension(RRspan(fs))>;
tof := func<v | &+[v[i]*monomialsY2[i] : i in [1..#monomialsY2]]>;

if IsEven(d) then
    harminv := ((d+3)*(d-1) + 3)/4;
else
    harminv := (d+3)*(d-1)/4;
end if;


//print "dim(H) = " cat " " cat IntegerToString(Integers() ! harminv);

print 2*harminv - d, harminv, harminv - d;

print RRdim(eqns), RRdim(eqnsp), RRdim(eqnsm);

tprime := ZeroMatrix(Rationals(),#monomialsY2,#monomialsY2);
for i,j,k,l in Integers(d) do  
    tprime[Index(monomialsY2,Y(i,j)*Y(k,l)),Index(monomialsY2,Y(i,l)*Y(k,j))] := 1;
end for;

prime := func<f|tof(RRv(f)*tprime)>;

// Finally seems to be working for general SICs. 
// Here's their dimensions ( Dim RR/(RRH + RRa) ) for small d:
// [d,dim] = [2,1], [3,5], [4,7], [5,11], [6, ], [7, ]

// Wow.  In d=2 (2 2 0), RRH = P1 meet P2, with
// P1 = (Y00 + Y01 + Y10 - Y11, Y01^2 + Y10^2)
// P2 = (Y00 - Y01 - Y10 - Y11, Y01^2 + Y10^2)

// Wow again.  In d=3 (3 3 0), RRH is prime.

// In d=4 (8 6 2) I ran out of memory computing the dimension.  May try another way. 

//S, toS := quo<RR|Idet>;

RH := ideal<R|[toR(f) : f in eqnsp]>;

RRa := ideal<RR|(d+1)*TrX - 1>;

RRY := ideal<RR|Y(0,0)-1>;

J := RRH + Idet + IT(2) + Ic;
Ja := J + RRa;
JY := J + RRY;
