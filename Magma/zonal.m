R<x,d> := PolynomialRing(Rationals(),2);
AssignNames(~R,["x","d"]);

function JacobiP(k,d,a,b,x)
    return HypergeometricSeries(-k,a+d-1,b,x);
end function;

function SQ(k)
    return GegenbauerPolynomial(k,d/2-1);
    //return ((d+2*k - 2)/(d-2))*Evaluate(GegenbauerPolynomial(k,(d-2)/2),x);
end function;
SQ(2);
