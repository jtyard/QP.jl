F2 := GF(2);
G := SL(2,F2);
q := Matrix(F2,[[1,1],[0,1]]);

#[g : g in G | Transpose(g)*q*g eq q];
// Why doesn't this give 6??  Should work

// On the other hand,

b := Matrix(F2,[[0,1],[1,0]]);
#[g : g in G | Transpose(g)*b*g eq b];

S3 := SymmetricGroup(3);
S2 := SymmetricGroup(2);
OW := WreathProduct(S3,S2);
CommutatorSubgroup(OW);

q00 := Matrix(F2,[[0,1],[0,0]]);
V00 := QuadraticSpace(q00);

q01 := Matrix(F2,[[0,1],[0,1]]);
V01 := QuadraticSpace(q01);

q10 := Matrix(F2,[[1,1],[0,0]]);
V10 := QuadraticSpace(q10);

q11 := Matrix(F2,[[1,1],[0,1]]);
V11 := QuadraticSpace(q11);

VH := QuadraticSpace(DirectSum(q11,q11));

OV00 := IsometryGroup(V00);
OV11 := IsometryGroup(V11);
OVH := IsometryGroup(VH);

SOVH := sub<OVH | [g : g in OVH | DicksonInvariant(VH,g) eq 0]>;
//CommutatorSubgroup(SOVH);

function OV(V)
    return IsometryGroup(V);
end function;

function SOV(V)
    O := OV(V);
    return sub<O|[g : g in O | DicksonInvariant(V,g) eq 0]>;
end function;   

function ThetaV(V)
    SO := SOV(V);
    return sub<SO|[g : g in SO | SpinorNorm(V,g) eq 0]>;
end function;  

function Comm(G)
    return CommutatorSubgroup(G);
end function;

function Orders(V)
    return [#Comm(SOV(V)), #Comm(OV(V)), #ThetaV(V), #SOV(V), #OV(V)];
end function;

print Orders(V00);
print Orders(V11);
print Orders(VH);

print Orders(QuadraticSpace(DirectSum(q00,q00)));

function qab(a,b)
    return Matrix(F2,[[a,1],[0,b]]);
end function;

function aut2(a,b,c,d)
    return IsometryGroup(QuadraticSpace(DirectSum(qab(a,b),qab(c,d))));
end function;