// This shows that the order of the EC stablizer is EXACTLY equal to the minimum of the orders of uf and -uf mod d.  It even gets the d=124 case right, which has an order-30 stabilizer (the first divisible by some prime other than 3 or 5)

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



for d:=4 to 48 do
    m, m0 := SquareFree((d-3)*(d+1));
    K := QuadraticField(m);
    OK := MaximalOrder(K);
    uf := FundamentalUnit(OK);
    if IsEven(d) then
	dd := 2*d;
    else
	dd := d;
    end if;
    if NormAbs(uf) eq -1 then
	up := uf^2;
    else
	if IsTotallyPositive(uf) then
	    up := uf;
	else
	    up := -uf;
	end if;
    end if;
    //print d, m, myorder(quo<OK|d*OK> ! up), myorder(quo<OK|d*OK> ! uf), myorder(quo<OK|dd*OK> ! up);
    ordp := Min(myorder(quo<OK|d*OK> ! up),myorder(quo<OK|d*OK> ! -up));
    ordf := Min(myorder(quo<OK|d*OK> ! uf),myorder(quo<OK|d*OK> ! -uf));
    ordp2 := myorder(quo<OK|d*OK> ! up);
    //if ord ne 3 and ord ne 6 then
    print d, m, ordp, ordp2, ordf;
    //end if;
end for;

// d:=5;
// %K := NumberField(QuadraticField((d-3)*(d+1)));
// %OK := RingOfIntegers(K);
