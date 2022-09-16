
// Assume d, r given, and that d is a prime 7 mod 12 with K having narrow class number 1.

Zd, ZZtoZd := Integers(d);

ad1 := PrimitiveElement(Zd);
a6r := ad1^(Integers() ! ((d-1)/(6*r)));
a3r := a6r^2;

print "Computing orbits";

orbs := {@ {@ (a3r^j)*i : j in [1..(3*r)]@} : i in [0..d-1]@};

N := #orbs;

mins := [Min(S) : S in orbs];
inds := [0 : i in [1..d]];
for i in [0..d-1] do
    inds[i+1] := Index([i in S : S in orbs],true);
end for;

names := ["z_" cat IntegerToString(Integers() ! i) : i in mins];

Q := Rationals();

print "creating " cat IntegerToString(N-1) cat "-dimensional projective space";

PZ := ProjectiveSpace(Q,N-1);
AssignNames(~PZ, names);
zz := [PZ.i : i in [1..N]];

z := func<i| zz[ inds[(Integers()!(Zd!i)) + 1] ] >;

norm2 := &+[z(i)^2 : i in [0..d-1]];

function delta(i)
    if (Integers(d) ! i) eq 0 then
	return 1;
    else
	return 0;
    end if;
end function;

ff := func<i,j|&+[z(k)*z(k+i)*z(k+j)*z(k+i+j) : k in [0..d-1]]  >;

f  := func<i,j| (d+1)*ff(i,j) - (delta(i) + delta(j))*norm2^2   >;

print "Defining equations";

ffeqns := [h : h in {@ ff(i,j) : j in [0..i], i in [0..(d+1)/2] @}];
eqns := [h : h in {@ f(i,j) : j in [0..i], i in [0..(d+1)/2] @}];

print "Making scheme";

proj := Scheme(PZ, eqns);
R := CoordinateRing(PZ);

print "Defining ideals";

RH := ideal<R|eqns>;
Ra := ideal<R|norm2 - 1>;
RHa := RH + Ra;


//time SS := Scheme(PZW, Radical(DefiningIdeal(S)));

function Del(f)
    return &+[Derivative(Derivative(f,z(k)),z(k)) : k in [0..d-1]];
end function;


//GL2 := GL(2,Integers(d));
//SL2 := SL(2,Integers(d));


//Qd := CyclotomicField(d);
//function beta(i,j)
//    return &+[Qd.1^(j*(a-b))*z(a)*z(a+i)*z(b)*z(b+i) : a,b in [0..d-1]];
//end function;

print "Defining tori and their ideals";

function  tt(a)
    t := ZeroMatrix(Q,N,N);
    for i in [1..N] do
	    t[i, inds[1+Integers() ! (a*(Zd!mins[i]))]] := 1;
	end for;
    return t;
end function;

function IT(a)
    eqnsa := {@ z(i) - z(i)^tt(a) : i in mins @};
    return ideal<R|eqnsa>;
end function;


function TT(a)
    return MatrixGroup<N,Q|tt(a)>;
end function;

monomials := MonomialsOfDegree(R,4);
N := #monomials;
R4 := VectorSpace(Q,N);

function Vf(f)
    return R4 ! [MonomialCoefficient(f,m) : m in monomials];
end function;

function Vspan(fs)
    return sub<R4|[Vf(f) : f in fs]>;
end function;

function Vdim(fs)
    return Dimension(Vspan(fs));
end function;

