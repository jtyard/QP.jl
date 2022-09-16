// d:=3; 

znames := ["z_" cat IntegerToString(i) : i in [0..d-1]];
wnames := ["w_" cat IntegerToString(i) : i in [0..d-1]];

names := znames cat wnames;
Q := Rationals();
// PZ<[z]> := ProjectiveSpace(Q,d-1);
// PW<[w]> := ProjectiveSpace(Q,d-1);
// PZW := DirectProduct(PZ,PW);

PZW := ProductProjectiveSpace(Q,[d-1,d-1]);
AssignNames(~PZW, names);
zz := [PZW.i : i in [1..d]];
ww := [PZW.i : i in [d+1..2*d]];

Zd, ZZtoZd := Integers(d);

z := func<i| zz[(Integers() ! (Zd ! i)) + 1]>;
w := func<i| ww[(Integers() ! (Zd ! i)) + 1]>;

R := CoordinateRing(PZW);

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

ffeqns := {@ ff(i,j) : i in [0..d-1], j in [0..d-1] | i ge j @};
eqns := {@ f(i,j) : i in [0..d-1], j in [0..d-1] | i ge j @};
eqns := [f : f in eqns];
S := Scheme(PZW, eqns);

//time SS := Scheme(PZW, Radical(DefiningIdeal(S)));

function Del(f)
    return &+[Derivative(Derivative(f,z(k)),w(k)) : k in [0..d-1]];
end function;

carpow := CartesianPower([0,1,2],d);
indices := [[A[k] : k in [1..d]]  : A in carpow | &+[A[k] : k in [1..d]] eq 2];
monomials := [&*[z(k)^A[k+1] : k in [0..d-1]] * &*[w(k)^B[k+1] : k in [0..d-1]]  : A,B in indices];

function fvec(i,j)
    fij := f(i,j);
    return [MonomialCoefficient(fij,m) : m in monomials];
end function;

function ffvec(i,j)
    ffij := ff(i,j);
    return [MonomialCoefficient(ffij,m) : m in monomials];
end function;

N := #monomials;

QN  := RModule(Integers(),N);
Mf  := sub<QN|[QN ! fvec(i,j)  : i,j in [0..d-1] | i le j]>;
Mff := sub<QN|[QN ! ffvec(i,j) : i,j in [0..d-1] | i le j]>;

if IsEven(d) then
    harminv := ((d+3)*(d-1) + 3)/4;
else
    harminv := (d+3)*(d-1)/4;
end if;

eqns2 := [ &+[(QN ! b)[i]*monomials[i] : i in [1..N]] : b in Basis(Mf)];

S2 := Scheme(PZW,eqns2);

print N, #eqns, #eqns2, harminv, Rank(Mf), Rank(Mff);

GL2 := GL(2,Integers(d));
SL2 := SL(2,Integers(d));

forbits := Matrix([[Index(eqns,f(i,j)) : j in [0..d-1]] : i in [0..d-1]]);

act := func<g,j| Eltseq(Vector(Integers(d), j)*g)>;


stab := [g : g in GL2 | Min([forbits[i+1,j+1] eq forbits[Integers() ! act(g,[i,j])[1]+1, Integers() ! act(g,[i,j])[2]+1] : i,j in [0..d-1] ]) eq true] ;

#stab;
//stab;

Qd := CyclotomicField(d);
function beta(i,j)
    return &+[Qd.1^(j*(a-b))*z(a)*w(a+i)*z(b)*w(b+i) : a,b in [0..d-1]];
end function;
