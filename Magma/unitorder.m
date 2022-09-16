function myorder(u)
    max_iters := 100;
    u0 := u^0;
    v := u^0;
    for k in [1..100] do
	v := v*u;
	if v eq u0 then
	    return k;
	end if;
    end for;
    "Order > 100.  Try increasing max_iters";
end function;

function unitorder(d)
    D := (d-3)*(d+1);
    D0 := FundamentalDiscriminant(D);
    K := QuadraticField(D0);
    OK := RingOfIntegers(K);
    uf := FundamentalUnit(OK);
    return myorder(quo<OK|5*OK> ! uf);
end function;

[unitorder(d) : d in [4..100]];
