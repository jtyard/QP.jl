// Stabilizer groupo
"loading";
S := [ScalarMatrix(IntegerRing(5), 2, 1),
      Matrix(IntegerRing(5), 2, 2, [ 0, 2, 2, 4 ]),
      Matrix(IntegerRing(5), 2, 2, [ 4, 3, 3, 0 ])];

// Coset representatives for the map s -> gl2/S
// Don't need anything else because cc acts trivially
// Therefore c acts as sig^4.
Gsigs := [ScalarMatrix(IntegerRing(5), 2, 1),
	  Matrix(IntegerRing(5), 2, 2, [ 3, 3, 3, 4 ]),
	  ScalarMatrix(IntegerRing(5), 2, 2),
	  Matrix(IntegerRing(5), 2, 2, [ 2, 3, 3, 3 ]),
	  ScalarMatrix(IntegerRing(5), 2, 4),
	  Matrix(IntegerRing(5), 2, 2, [ 2, 2, 2, 1 ]),
	  ScalarMatrix(IntegerRing(5), 2, 3),
	  Matrix(IntegerRing(5), 2, 2, [ 3, 4, 4, 1 ])];

h_of_g := func<g | Matrix(IntegerRing(5), 2, 2, [ 1, 0, 0, Determinant(g) ])>;

Hsigs := [h_of_g(g) : g in Ssigs];
Ssigs := [Gsigs[i]*Hsigs[i]^-1 : i in [1..#Ssigs]];
"loaded";


// Test cocycle condition 
//for g0,g1 in S do
//    g0,g1;
//end for;
