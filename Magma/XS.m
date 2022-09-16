F := CyclotomicField(8);
w := F.1;
i := w^2;
s2 := w + w^-1;

I := Matrix(F,2,[1,0,0,1]);
X := Matrix(F,2,[0,1,1,0]);
Y := Matrix(F,2,[0,-i,i,0]);
S := Matrix(F,2,[1,0,0,i]);
H := Matrix(F,2,[1,1,1,-1])/(w+w^-1);
Z := S^2;


D4 := MatrixGroup<2,F|X,Z>;
Q8 := MatrixGroup<2,F|i*X,i*Z>;
XS := MatrixGroup<2,F|i*X,S/w>;
P2 := MatrixGroup<2,F|i*I,X,Z>;
#XS;



