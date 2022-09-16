for d in [4..20] do
    D := (d-3)*(d+1);
    K := QuadraticField(D);
    OK := MaximalOrder(K);
    u := (d-1 + Sqrt(K ! D))/2;
    OD := sub<OK | u>;
    print d, FundamentalUnit(OD) eq u;
    //if Norm(FundamentalUnit(OK)) eq -1 and Norm(FundamentalUnit(OD)) eq -1 then
    //    print d, D;
    //end if;
end for;
