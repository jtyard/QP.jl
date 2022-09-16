dims := [d : d in PrimesUpTo(50) | IsOdd(d) and d mod 3 eq 1];
//dims := [5779];
print dims;

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

function fundamental_unit(K)
    uf := FundamentalUnit(K);
    u1 := Evaluate(uf,RealPlaces(K)[1]);

    if u1 lt 0 then
        uf := -uf;
    end if;

    if u1^2 lt 1 then
        uf := 1/uf;
    end if;

    return uf;
end function;   
  


for d in dims do
    D := (d-3)*(d+1);
    K := QuadraticField(D);
    OK := MaximalOrder(K);
    //u := (d-1 + Sqrt(K ! D))/2;
    uf := fundamental_unit(K);
    
    h := ClassNumber(OK);
    
    if Norm(uf) eq 1 then
        up := uf;
        hp := 2*h;
    else
        hp := h;
        up := uf^2;
    end if;
    
    m,f := Squarefree(D);

    fac := Factorization(ideal<OK|d>);
    P1 := fac[1][1];
    P2 := fac[2][1];
    
    n1 := [#HeckeCharacterGroup(P1), #HeckeCharacterGroup(P1,[1]), #HeckeCharacterGroup(P1,[2]), #HeckeCharacterGroup(P1,[1,2])];
    n2 := [#HeckeCharacterGroup(P2), #HeckeCharacterGroup(P2,[1]), #HeckeCharacterGroup(P2,[2]), #HeckeCharacterGroup(P2,[1,2])];
    n := [#HeckeCharacterGroup(ideal<OK|d>), #HeckeCharacterGroup(ideal<OK|d>,[1]), #HeckeCharacterGroup(ideal<OK|d>,[2]), #HeckeCharacterGroup(ideal<OK|d>,[1,2])];

    r1 := myorder(quo<OK|P1> ! uf);
    r2 := myorder(quo<OK|P2> ! uf);
    r12  := myorder(quo<OK|d> ! uf);

    r := myorder(quo<OK|d> ! up)/3;

    e1 := myorder(quo<OK|P1> ! -1);
    e2 := myorder(quo<OK|P2> ! -1);
    e  := myorder(quo<OK|d> ! -1);



    if true then
        print d,m,h,hp,Norm(uf);
        print d mod 4;
        print n1;
        print n2;
        print n;
        print h*2*(d-1)^2/r12, (d-1)^2*hp/(3*r),(d-1)*hp/(3*r);
        print [r1,r2,r12,r];
        print [e1,e2,e];
        print "";
        //OD := sub<OK | u>;
        //print d, FundamentalUnit(OD) eq u;
    end if;

    print "* ", n[4], n1[4]*(d-1), n2[4]*(d-1),hp*(d-1)^2/(3*r);

end for;
