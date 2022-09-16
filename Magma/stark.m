d:=5;
precN := 70;


if IsEven(d) then
    dd := 2*d;
else
    dd := d;
end if;

D := (d-3)*(d+1);

K<s> := QuadraticField(D);
OK := MaximalOrder(K);
KX<x> := PolynomialRing(K);
inf := InfinitePlaces(K);


Chi1 := HeckeCharacterGroup(ideal<OK|dd>, [1]);
G1, ph1 := RayClassGroup(ideal<OK|dd>, [1]);
deg := #G1;
sig1 := ph1(G1.1);
Chis1 := [Chi1.1^i : i in [0..deg-1]];

Chi12 := HeckeCharacterGroup(ideal<OK|dd>, [1,2]);
G12, ph12 := RayClassGroup(ideal<OK|dd>, [1,2]);

sig12 := Chi12.1;
ep12 := Chi12.2;
Chis12 := [sig12^i*ep12^j : i in [0..deg-1],j in [0..1]];

//c := ph1(4*sig1);


//chi1 := Chi.1;
//chi2 := Chi.2;
//Chi12 := [chi1^i*chi2^j : i in [0..7], j in [0,1]];


print [Norm(Conductor(chi)) : chi in Chis1];
print [chi(sig1^4) : chi in Chis1];  
print [IsPrimitive(chi) : chi in Chis1];


function Lval(chi)
    if chi(sig1^(Integers()! (deg/2))) eq 1 then
        return 0;
    end if;
    return Evaluate(LSeries(chi : Precision:=precN),0 : Derivative:=1);
end function;

// Uncomment this to recompute the L-values
Lvals1 := [Lval(chi) : chi in Chis1];

zetavals1 := [&+[( ComplexField(precN) ! (Chis1[k]^-1)(ph1(g)))*Lvals1[k] : k in [1..deg]]/#G1 : g in G1];

// Recognize element of the quadratic field K
function recognize(a)
    f := PowerRelation(a,2);
    fac := Factorization(PolynomialRing(K) ! f);
    vals := [Abs(Evaluate(PolynomialRing(RealField(Precision(a))) ! ff[1],a)) : ff in fac];
    minval, minj := Min(vals);
    return -Coefficient(fac[minj][1],0);
end function;

function topoly(galorbit)
    y := PolynomialRing(RealField(precN)).1;
    f := &*[y-a : a in galorbit]; 
    return KX ! [recognize(b) : b in Coefficients(f)];
end function;

starks1 := [Exp(2*Real(a)) : a in zetavals1];
starks2 := [Exp(Real(a)) : a in zetavals1];
starks3 := [Exp(Real(a))*Sqrt(RealField(precN) ! (d+1)) : a in zetavals1];

f:=topoly(starks1); print f;
print Factorization(Evaluate(f,x^2/(d+1)));
//print topoly(starks2);
//print topoly(starks3);

//Ff := NumberField(f);
//F12 := NumberField(RayClassField(Chi12));
//F1 := NumberField(RayClassField(Chi1));
//F2 := NumberField(RayClassField(HeckeCharacterGroup(ideal<OK|dd>, [2])));



//[Discriminant(NumberField(Polynomial(IntegerRelation([a^i : i in [0,1,2]])))) : a in A];

