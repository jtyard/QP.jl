// rcfactor(-1,5,[1,2]) e.g.

function classnumber(d)
    return ClassNumber(QuadraticField(d-3)*(d+1));
end function;



// Factor d in the rcf of Q((d-3)(d+1)) with conductor m inf, inf subset [1,2]
function factord(d,m,inf)
    K := QuadraticField((d-3)*(d+1));
    OK := RingOfIntegers(K);
    //"Creating ray class group";
    rayclassgroup, representing_ideal := RayClassGroup(m*OK,inf);

    //"Building abelian extension";
    FA := AbelianExtension(representing_ideal);
    
    //"Constructing number field";
    F:= NumberField(FA);
    OF:= RingOfIntegers(F);
    return Factorization(d*OF);
end function;

function efg(d,m,inf)
    //"help me";
    fac := factord(d,m,inf);
    q := Norm(AbsoluteNorm(d*RingOfIntegers(QuadraticField((d-3)*(d+1)))));
    return [fac[1][2],Round(Log(q,AbsoluteNorm(fac[1][1]))),#fac];
end function;

function dd(d)
    if IsEven(d) then
	return 2*d;
    else
	return d;
    end if;
end function;
	  
	 
