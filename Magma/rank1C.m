

names := ["Y" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];

Q := Rationals();

Cd := CyclotomicField(2*d);
zdd := Cd.1;
zd := zdd^2;

PY := ProjectiveSpace(Cd,d^2-1);
AssignNames(~PY, names);
YY := [PY.i : i in [1..d^2]];

Zd, ZZtoZd := Integers(d);

Y := func<i,j| YY[d*(Integers() ! (Zd ! i)) + (Integers() ! (Zd ! j)) + 1]>;

Ymat := Matrix([[Y(i,j) : j in [0..d-1]] : i in [0..d-1]]); Y2mat := KroneckerProduct(Ymat,Ymat);

Tr1 := Trace(Ymat);
Tr2 := Trace(Ymat^2);

Y0 := func<i,j|Y(i,j)-Tr1>;

norm2 := &+[Y(i,i) : i in [0..d-1]];

Del := func<f|&+[Derivative(f,Y(a,a)) : a in [0..d-1]]>;


function delta(i)
    if (Integers(d) ! i) eq 0 then
	return 1;
    else
	return 0;
    end if;
end function;

function delta2(a,b)
    if (Integers(d) ! a) eq 0 and (Integers(d) ! b) eq 0 then
	return 1;
    else
	return 0;
    end if;
end function;

// Bad notation since these are the FTs of the Bj. 
fj := func<i,j|&+[Y(a,a+i)*Y(a+i+j,a+j) : a in [0..d-1]]>;

//Cj := func<i,j|&+[Y(a,a+i)*Y(a+i+j,a+j) : a in [0..d-1]] >;

if IsEven(d) then
    h20p := (d/2)^2 + d/2;
    h20m := (d/2)^2 - d/2; 
else
    h20p := (d+3)*(d-1)/4;
    h20m := (d-3)*(d+1)/4;
end if;

b0p := (Tr1 + Tr2)/2;
b0m := (Tr1 - Tr2)/2;

//Bp := {@ &+[Y0(a,a+j1)*Y0(a+j1+j2,a+j2)+Y0(a,a+j2)*Y0(a+j1+j2,a+j1) : a in [0..d-1]] : j1,j2 in [0..d-1]@};
//Bm := {@ &+[Y0(a,a+j1)*Y0(a+j1+j2,a+j2)-Y0(a,a+j2)*Y0(a+j1+j2,a+j1) : a in [0..d-1]] : j1,j2 in [0..d-1]@};
//
//Bp := [f : f in Bm | f ne 0];
//Bm := [f : f in Bm | f ne 0];


fjp := func<a,b|(fj(a,b)+fj(b,a))/2>;
Hp := func<i,j|fjp(i,j)>;	

fjm := func<a,b|(fj(a,b)-fj(b,a))/2>;
// keep guessing here... 
Hm := func<i,j|fjm(i,j)>;


Aj := func<j1,j2 | (-zdd)^(j1*j2)*(&+[zd^(j2*a)*Y(a,a+j1) : a in [0..d-1]])>;
Bj := func<j1,j2| Aj(j1,j2)*Aj(-j1,-j2)>;

Bftj := func<j1,j2|&+[Bj(j1,a)*zd^(a*j2) : a in [0..d-1]]/d>;
Cj := func<j1,j2|&+[Y(a,a+j1)*Y(a+j1+j2,a+j2) : a in [0..d-1]]>;

Cjh := func<j1,j2|Cj(j1,j2) - ((d*delta(j1) - delta(j2))*Tr1^2 + (-delta(j1) + d*delta(j2))*Tr2)/(d^2-1)>;

Dj := func<j1,j2|&+[Cj(a,j2)*zd^(-a*j1) : a in [0..d-1]]>;


// Define generalized Pauli matrices
Xd:=MatrixRing(Cd,d)!0;
Xd[1,d]:=1;
for i:=1 to d-1 do
    Xd[i+1,i]:=1;
end for;

Zd := DiagonalMatrix(Cd,[zd^i:i in [0..d-1]]);

XZ := func<j| Xd^(Integers() ! j[1])*Zd^(Integers() ! j[2]) >;

Delta := func<j| (-zdd)^(Integers() ! (j[1]*j[2]))*XZ(j) >;

Ymat := Matrix([[Y(i,j) : j in [0..d-1]] : i in [0..d-1]]);
Y2mat := KroneckerProduct(Ymat,Ymat);

Idmat := MatrixRing(Q,d) ! 1;
va := func<a| Transpose(Matrix(Idmat[(a mod d)+1]))>;

SWAP := &+[KroneckerProduct(va(a)*Transpose(va(b)),va(b)*Transpose(va(a))):a,b in [0..d-1]];

Dj := func<j1,j2|Trace(Y2mat*SWAP*KroneckerProduct(Delta([j1,j2]),Delta([-j1,-j2])))>;

Bjh := func<j1,j2| (1-delta2(j1,j2))*(Bj(j1,j2) + (Tr1^2 - d*Tr2)/(d^2-1))>;

Bjp := func<j1,j2| (Bj(j1,j2) + Dj(j1,j2))/2>; 
Bjph := func<j1,j2| Bjp(j1,j2) - (delta2(j1,j2)*d + 1)*Bjp(0,0)/(d+1)>;

Bjm := func<j1,j2| (Bj(j1,j2) - Dj(j1,j2))/2>;
Bjmh := func<j1,j2| Bjm(j1,j2) - (delta2(j1,j2)*d - 1)*Bjm(0,0)/(d-1)>;

// Define the indexing group and module
sl2 := SL(2,Integers(d));
gl2 := GL(2,Integers(d));
M := GModule(gl2);
w := gl2 ! [0,1,1,0];

ind2 := [[Integers() ! jj[1][1], Integers() ! jj[1][2]] : jj in {@ {@ j, -j, j*w, -j*w @} : j in M | j ne 0 @}];


lightcone := [Cjh(a,a) : a in [0..(d-1)/2]] cat [Cjh(a,-a) :a in [1..(d-1)/2]];
eqnsp := [Cjh(j1,j2) + Cjh(j2,j1) : j1 in [0..(d-1)/2], j2 in [-(d-1)/2..(d-1)/2] | Abs(j2) le j1];
eqnsm := [Cjh(j1,j2) - Cjh(j2,j1) : j1 in [1..(d-1)/2], j2 in  [-(d-1)/2..(d-1)/2] | Abs(j2) lt j1];

//eqnsp := {@  Cjh(j[1],j[2]) + Cjh(j[2],j[1]) : j in ind2 @};
//eqnsm := {@  Cjh(j[1],j[2]) - Cjh(j[2],j[1]) : j in ind2 | j[1] ne j[2] @};

"#eqnsp, #eqnsm=", #eqnsp, #eqnsm;
"h20p, h20m=", h20p, h20m;
//"[h20p, h20m,h20p+h20m]=" , [h20p, h20m,h20p+h20m];

R := CoordinateRing(PY);
minors := [Y(j1,k1)*Y(j2,k2) - Y(j1,k2)*Y(j2,k1) : j1,j2,k1,k2 in [0..d-1] | j1 lt j2 and k1 lt k2];
I2 := ideal<R|minors>;
//2*d-1;

Ip := ideal<R|eqnsp>;
Im := ideal<R|eqnsm>;

IH := ideal<R|eqnsp,eqnsm>;
//Dimension(IH);

IH2 := ideal<R|eqnsp,eqnsm,minors>;
Icone := ideal<R|lightcone>;
Rp := ideal<R|YY>;

Bmat := func<j1,j2|(KroneckerProduct(Delta([j1,j2]),Delta([-j1,-j2])) + KroneckerProduct(Delta([-j1,-j2]),Delta([j1,j2])))>;

Dmatm := func<j1,j2|SWAP*(KroneckerProduct(Delta([j1,j2]),Delta([-j1,-j2])) + KroneckerProduct(Delta([-j1,-j2]),Delta([j1,j2])))>;
