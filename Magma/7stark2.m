//d := 7;

dims := [ 7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109, 127, 139, 151, 157, 163, \
    181, 193, 199, 211, 223, 229, 241, 271, 277, 283, 307, 313, 331, 337, 349, 367,  \
    373, 379, 397, 409, 421, 433, 439, 457, 463, 487, 499, 523, 541, 547, 571, 577, \
    601, 607, 613, 619, 631, 643, 661, 673, 691, 709, 727, 733, 739, 751, 757, 769, \
    787, 811, 823, 829, 853, 859, 877, 883, 907, 919, 937, 967, 991, 997, 1009, \
    1021, 1033, 1039, 1051, 1063, 1069, 1087, 1093, 1117, 1123, 1129, 1153, 1171,  \
    1201, 1213, 1231, 1237, 1249, 1279, 1291, 1297, 1303, 1321, 1327, 1381, 1399,  \
    1423, 1429, 1447, 1453, 1459, 1471, 1483, 1489, 1531, 1543, 1549, 1567, 1579,  \
    1597, 1609, 1621, 1627, 1657, 1663, 1669, 1693, 1699, 1723, 1741, 1747, 1753,  \
    1759, 1777, 1783, 1789, 1801, 1831, 1861, 1867, 1873, 1879, 1933, 1951, 1987,  \
    1993, 1999 ];

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




//uf := 1+s;
//up := uf^2;

P := ideal<OK|d>;
fac := Factorization(P);
P1 := fac[1][1];
P2 := fac[2][1];
 
Chi0 := HeckeCharacterGroup(ideal<OK|1>,[1,2]);
gen0 := [chi : chi in Generators(Chi0)];
G0 := RayClassGroup(ideal<OK|1>,[1,2]);


Chi1 := HeckeCharacterGroup(P1,[1,2]); 
gen1 := [chi : chi in Generators(Chi1)];
G1 := RayClassGroup(P1,[1,2]);

Chi2 := HeckeCharacterGroup(P2,[1,2]); 
gen2 := [chi : chi in Generators(Chi2)];
G2 := RayClassGroup(P2,[1,2]);

Chi12 := HeckeCharacterGroup(P2,[1,2]); 
G12 := RayClassGroup(ideal<OK|d>,[1,2]);

hp := #Chi0;
print #Chi0, #Chi1, #Chi2;

Kinf := NumberField(RayClassField(Chi0));
OKinf := MaximalOrder(Kinf);

print "h, hp = ", ClassNumber(K), hp;

print "Chi0: ", [Order(chi) : chi in Generators(Chi0)]; 
for chi in Generators(Chi0) do
    print Conductor(chi);
end for;

print "Chi1:", [Order(chi) : chi in Generators(Chi1)]; 

for chi in Generators(Chi1) do
    print Conductor(chi);
end for;

print "Chi2:", [Order(chi) : chi in Generators(Chi2)]; 

for chi in Generators(Chi2) do
    print Conductor(chi);
end for;

// Recognize element of the quadratic field K (later fix to do Kinf..)
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

function Lval(chi)
    return Evaluate(LSeries(AssociatedPrimitiveCharacter(chi) : Precision:=precN),0 : Derivative:=1);
end function;


//a := Evaluate(LSeries(Chi1.1 : Precision:=precN),0 : Derivative:=1);

//print topoly([Exp(a),Exp(-a)]);

//print topoly([Exp(2*a),Exp(-2*a)]);

// Holy shit the negative of my magic unit is a stark unit (duh it should have been positive!) 
// min poly x^2 - (s + 1)x + 1
