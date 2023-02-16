d := 7;
precN := 70;
RealN := RealField(precN);

function myorder(u)
    nmax := #MultiplicativeGroup(Parent(u));
    u0 := u^0;
    ui := u;
    for i:=1 to nmax do
	if ui eq u0 then
	    return i;
	end if;
	ui := ui * u;
    end for;
    return 0;
end function;

D := (d-3)*(d+1);
K<s> := QuadraticField(D);
OK := MaximalOrder(K);
KX<x> := PolynomialRing(K);


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

uf := 1+s;
up := uf^2;

P1 := ideal<OK|3-s>;
P2 := ideal<OK|3+s>;
 
Chi := HeckeCharacterGroup(P1,[1]);
G := RayClassGroup(P1,[1]);

a := Evaluate(LSeries(Chi.1 : Precision:=precN),0 : Derivative:=1);

print topoly([Exp(a),Exp(-a)]);

print topoly([Exp(2*a),Exp(-2*a)]);

// The negative of my magic unit is a stark unit (duh it should have been positive!) 
// min poly x^2 - (s + 1)x + 1
