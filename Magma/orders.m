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
D0 := FundamentalDiscriminant(D);

f := Integers() ! Sqrt(D/D0);

K := QuadraticField(D);
OK := MaximalOrder(K);

u := (d-1 + Sqrt(K ! D))/2;

up := FundamentalUnit(K);
if Norm(up) eq -1 then
    up := up^2;
end if;

rmax := Integers() ! (myorder(quo<OK|d*OK> ! up)/3);

print rmax;

orders := [sub<OK | ff*OK.2> : ff in Divisors(f)];

print Divisors(f);
print [rmax/Min([r : r in [1..rmax] | up^r in O]) : O in orders];
print [myorder(quo<OK|d*ff*OK> ! up)/3 : ff in Divisors(f)];

print 1 + up^rmax + up^-rmax;

OD := sub<OK | u>;

//if myorder(quo<OD|d*OD> ! u) ne 3 then
//    print d, D,  myorder(quo<OD|d*OD> ! u);
//end if;
    //if Norm(FundamentalUnit(OK)) eq -1 and Norm(FundamentalUnit(OD)) eq -1 then
    //    print d, D;
    //end if;
