// Shit yeah.  This shows that the order of the EC stablizer is EXACTLY equal to the minimum of the orders of uf and -uf mod d.  
// It even gets the d=124 case right, which has an order-30 stabilizer (the first divisible by some prime other than 3 or 5)
// Also looks like the order mod dd is *always* twice the order mod d (when d is even). Note that this the orders can change mod 2d if d is odd.

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



for d:=4 to 124 do
    m, m0 := SquareFree((d-3)*(d+1));

    if IsEven(d) then
	dd := 2*d;
    else
	dd := d;
    end if;

    K := QuadraticField(m);
    OK := MaximalOrder(K);
   
    uf := FundamentalUnit(OK);

    if Real(uf) le 0 then
	uf := -uf;
    end if;
    
    if NormAbs(uf) eq -1 then
	up := uf^2;
    else
	up := uf;
    end if;
    //print d, m, myorder(quo<OK|d*OK> ! up), myorder(quo<OK|d*OK> ! uf), myorder(quo<OK|dd*OK> ! up);

    ordufd := myorder(quo<OK|d*OK> ! uf);
    ordufdd := myorder(quo<OK|dd*OK> ! uf^2);
    ordupd := myorder(quo<OK|d*OK> ! up);
    ordupdd := myorder(quo<OK|dd*OK> ! up^2);
    
    if true then
	print d, NormAbs(uf), "|",  ordufd, ordupd, "|",  ordufdd, ordupdd;
    end if;
   
end for;

// d:=5;
// %K := NumberField(QuadraticField((d-3)*(d+1)));
// %OK := RingOfIntegers(K);
