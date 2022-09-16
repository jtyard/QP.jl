for m in M do
    m;
    phi2 := XZ(m)*phi;
    ff:= map<M -> F | m :-> ip(phi2,XZ(m)*phi2)/phinorm >;
    ff2:= map<M -> F | m :-> ip(phi2,Delta(m)*phi2)/phinorm >;
    #{ff(m) : m in M};
    #{ff2(m) : m in M};
    //SS := [g : g in sl2 | equal(act(g)*ff,ff)];
    //SS2 := [g : g in sl2 | equal(act(g)*ff2,ff2)];
    
    //if #SS gt 1 then
	//centerS := SS;
	//centerm := m;
    //end if;

    //if #SS2 gt 1 then
	//centerS2 := SS2;
        //centerm2 := m;
    //end if;

end for;
