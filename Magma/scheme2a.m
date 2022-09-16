d:=2;
load "scheme2.m";

F := CyclotomicField(12);
C4 := CyclotomicField(4);
C3 := CyclotomicField(3);

CS := S;

rationalpoints := RationalPoints(CS,F);
points := [Eltseq(p) : p in rationalpoints];

PS := IrreducibleComponents(S);
PS3 := IrreducibleComponents(BaseChange(S,C3));
PS4 := IrreducibleComponents(BaseChange(S,C4));

ps :=  [RationalPoints(X,F) : X in PS];
ps4 := [RationalPoints(X,F) : X in PS4];
ps3 := [RationalPoints(X,F) : X in PS3];

print [[Index(points,Eltseq(p)) : p in A ] : A in ps];

//print ps4;
//print ps3;

sic1 := [points[1], points[2], points[5], points[6]];
sic2 := [points[3], points[4], points[7], points[8]];

prod := func<p,q|Trace(FtoM(p)*FtoM(q))>;
print Matrix([[prod(p,q) : q in sic1] : p in sic1]);
print Matrix([[prod(p,q) : q in sic2] : p in sic2]);
print [prod(A[1],A[2]) : A in ps4];
print [prod(A[1],A[2]) : A in ps3];

sigX := Matrix([[0,1],[1,0]]);
sigZ := Matrix([[1,0],[0,-1]]);
ii := C4.1;
sigY := Matrix([[0, -ii],[ii,0]]);
FtoV := func<p |  Matrix(F, FtoM(p))>;
print [[Trace(sigX*FtoV(p)), Trace(sigY*FtoV(p)), Trace(sigZ*FtoV(p))] : p in points];


sl2 := SL(2,GF(2));
M := GModule(sl2);

v00 := M ! [0,0];
v01 := M ! [0,1];
v10 := M ! [1,0];
v11 := M ! [1,1];

P0 := FtoM(sic1[1]);
ZZ := Integers();

P := func<j | sigX^(ZZ!j[1])*sigZ^(ZZ!j[2])*P0*sigZ^(ZZ!j[2])*sigX^(ZZ!j[1])>;

Tr3 := func<i,j,k|Trace(P(i)*P(j)*P(k))>;
Tr2 := func<i,j|Trace(P(i)*P(j))>;
