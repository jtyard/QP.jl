names := ["Y" cat IntegerToString(i) cat IntegerToString(j) : i,j in [0..d-1]];

Q := Rationals();

//Cd := CyclotomicField(2*d);
//zdd := Cd.1;
//zd := zdd^2;

PY := ProjectiveSpace(Q,d^2-1);
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

Cj := func<j1,j2|&+[Y(a,a+j1)*Y(a+j1+j2,a+j2) : a in [0..d-1]]>;

Cjh := func<j1,j2|Cj(j1,j2) - ((d*delta(j1) - delta(j2))*Tr1^2 + (-delta(j1) + d*delta(j2))*Tr2)/(d^2-1)>;


Ymat := Matrix([[Y(i,j) : j in [0..d-1]] : i in [0..d-1]]);
Y2mat := KroneckerProduct(Ymat,Ymat);

Idmat := MatrixRing(Q,d) ! 1;
va := func<a| Transpose(Matrix(Idmat[(a mod d)+1]))>;

SWAP := &+[KroneckerProduct(va(a)*Transpose(va(b)),va(b)*Transpose(va(a))):a,b in [0..d-1]];


sl2 := SL(2,Integers(d));
gl2 := GL(2,Integers(d));
M := GModule(gl2);
w := gl2 ! [0,1,1,0];

ind2 := [[Integers() ! jj[1][1], Integers() ! jj[1][2]] : jj in {@ {@ j, -j, j*w, -j*w @} : j in M | j ne 0 @}];


lightcone := [Cjh(a,a) : a in [0..(d-1)/2]] cat [Cjh(a,-a) :a in [1..(d-1)/2]];
forw :=  [Cjh(j1,j2) + Cjh(j2,j1) : j1 in [1..(d-1)/2], j2 in  [-(d-1)/2..(d-1)/2] | Abs(j2) lt j1];	  
eqnsp := [Cjh(j1,j2) + Cjh(j2,j1) : j1 in [0..(d-1)/2], j2 in [-(d-1)/2..(d-1)/2] | Abs(j2) le j1];
eqnsm := [Cjh(j1,j2) - Cjh(j2,j1) : j1 in [1..(d-1)/2], j2 in  [-(d-1)/2..(d-1)/2] | Abs(j2) lt j1];

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
IH2 := IH + ideal<R|Tr1^2 - Tr2>;
IH3 := ideal<R|eqnsp,minors>;


Icone := ideal<R|lightcone>;
Iforw := ideal<R|forw>;
// The irrelevant ideal
Rp := ideal<R|YY>;

Bmat := func<j1,j2|(KroneckerProduct(Delta([j1,j2]),Delta([-j1,-j2])) + KroneckerProduct(Delta([-j1,-j2]),Delta([j1,j2])))>;
Dmatm := func<j1,j2|SWAP*(KroneckerProduct(Delta([j1,j2]),Delta([-j1,-j2])) + KroneckerProduct(Delta([-j1,-j2]),Delta([j1,j2])))>;
