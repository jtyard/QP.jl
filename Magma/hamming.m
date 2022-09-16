n:= 4;
Q := Rationals();
A := AffineSpace(Q,n);
R := CoordinateRing(A);

gens := [R.i^1-1 : i in [1..n]];
I := ideal<R|gens>;

// Oops the dimension is zero, how to find number of generators?