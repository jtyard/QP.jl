function myorder(u)
    nmax := 100;
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
K := QuadraticField(D);
OK := MaximalOrder(K);
D0 := FundamentalDiscriminant(D);
f := Integers() ! Sqrt(D/D0);

divs := Divisors(f);

uf := FundamentalUnit(OK);
uf1 := Evaluate(uf,RealPlaces(K)[1]);
if 0 lt uf1 and uf1 lt 1 then
    uf := 1/uf;
else 
    if -1 lt uf1  and uf1 lt 0 then
        uf := -1/uf;
    else 
        if uf1 lt -1 then
            uf := -uf;
        end if;
    end if;
end if;

if Norm(uf) eq 1 then
    up := uf;
else 
    up := uf^2;
end if;

r := Integers() ! (myorder(quo<OK|d*OK> ! up)/3);

u3 := (d-1 + Sqrt(K ! D))/2;
OD := sub<OK | u3>;

print d, Norm(uf), r,D0, D, Discriminant(OD);
print divs;
print Divisors(r);

x := PolynomialRing(K).1; 

if not IsSquare(up) then
    E := ext<K|x^2 - up>;
    print ideal<OK|Discriminant(E)> eq ideal<OK|4>;
else 
    print "Sqrt(up) already contained in K";
end if; 

if not IsSquare(-up) then
    E2 := ext<K|x^2 + up>;
    print Discriminant(E2), Norm(Discriminant(E2));
    print ideal<OK|Discriminant(E2)> eq ideal<OK|4>;
else 
    print "Sqrt(-up) already contained in K";
end if; 
