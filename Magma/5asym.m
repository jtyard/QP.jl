
//Ssig := [g : g in sl2 | equal(act(g)*f*sig,f)];
//Ssig2 := [g : g in sl2 | equal(act(g)*f2*sig,f2)];

//Sc := [g : g in sl2 | equal(act(g)*f*c,f)];
//Sc2 := [g : g in sl2 | equal(act(g)*f2*c,f2)];

//S0 := [g : g in sl2 | equal(act(g)*f,f)];
//S02 := [g : g in sl2 | equal(act(g)*f2,f2)];

j := Random(M);
k := Random(M);
j,k;


cocycle:=func<j,k | (-z2d)^(Integers() ! (j[2]*k[1] - j[1]*k[2]))>;
pairing:=func<j,k | zd^(Integers() ! (j[2]*k[1] - j[1]*k[2]))>;


pairing(j,k) eq cocycle(j,k)/cocycle(k,j);
Delta(j)*Delta(k) eq cocycle(j,k)*Delta(j+k);
Delta(j)*Delta(k)*Delta(-j)*Delta(-k) eq pairing(j,k)*Delta(M!0);


function sum(S)
    s := 0;
    for a in S do
	s := s+a;
    end for;
    return s;
end function;

function prod(S)
    p := 1;
    for a in S do
	p := p*a;
    end for;
    return p;
end function;


//P := sum({f2(j)*Delta(-j)/d : j in M});

stab:=func<f | [g : g in sl2 | equal(act(g)*f,f)]>;
// Now this works, giving 3 elements
// stab(f2);
//sl2stab := [<g,s> : g in sl2, s in Gal | equal(act(g)*f2*s,f2)];
//"#sl2stab ", #sl2stab;

"Computing sigstab";
sigstab := [<g,i> : g in gl2, i in [1..8] | equal(act(g)*f2*(sig^i),f2)];

//#{(sig^j)(f2([1,1])) : j in [0..7]};
//#{f2(j) : j in M};
//{(sig^j)(f2([1,1])) : j in [0..7]} subset {f2(j) : j in M};


//for g in sl2 do
//    for h in sl2 do
//	if not equal(act(g*h)*f2,act(g)*(act(h)*f2)) then
//	    g,h;
//	end if;
//    end for;
//end for;
    
    
	
	

	
		   

