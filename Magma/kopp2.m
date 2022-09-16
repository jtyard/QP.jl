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


for d in [4..10000] do
    D := (d-3)*(d+1);
    K := QuadraticField(D);
    OK := MaximalOrder(K);
    u := (d-1 + Sqrt(K ! D))/2;
    OD := sub<OK | u>;
    if myorder(quo<OD|d*OD> ! u) ne 3 then
	print d, D,  myorder(quo<OD|d*OD> ! u);
    end if;
    //if Norm(FundamentalUnit(OK)) eq -1 and Norm(FundamentalUnit(OD)) eq -1 then
    //    print d, D;
    //end if;
end for;
