


RealCyclotomicField := function(n)
    F := CyclotomicField(n);
    return sub<F|F.1 + F.1^-1>;
end function;

Dn := function(n)
    if (n mod 4) eq 0 then 
        return 1;
    else
        return 4 - (RealCyclotomicField(n).1)^2;
    end if;
end function;

su2levelk := function(k) 
    K := RealCyclotomicField(k+2);
    return QuaternionAlgebra(K,-1 - K.1,Dn(k+2));
end function;


A := su2levelk(3);

OA := MaximalOrder(A);

allinvs := function(OO)    
    return [#TwoSidedIdealClassGroup(OO),#ConjugacyClasses(OO), #LeftIdealClasses(OO),#RightIdealClasses(OO)];
end function;