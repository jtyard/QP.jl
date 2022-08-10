// n := 3 // gives Clifford group


SetClassGroupBounds("GRH");


F := CyclotomicField(2^n);
K := sub<F|F.1 + F.1^-1>;
A := QuaternionAlgebra(K,-1,-1);
OA := MaximalOrder(A);
OL := Order(Basis(A));
ii := Basis(A)[2];
jj := Basis(A)[3];
kk := Basis(A)[4];
s := Sqrt(K ! 2);
OH := Order([1,(ii + 1)/s ,(jj + 1)/s, (1 + ii + jj + kk)/2]);

// Or create a non-maximal one??
// OK := MaximalOrder(K);
// OL := Order(OK,Basis(A));

invs := function(OO)    
    return [#TwoSidedIdealClassGroup(OO),#ConjugacyClasses(OO), Mass(OO)];
end function;

allinvs := function(OO)    
    return [#TwoSidedIdealClassGroup(OO),#ConjugacyClasses(OO), #LeftIdealClasses(OO),#RightIdealClasses(OO),Mass(OO)];
end function;

// Strange that Left/Right classes don't work for for n geq 5 because 
// because right_ideal_classes = 2_sided_ideal_classes x conjugacy_classes
// i.e. conjugacy_classes = right_ideal_classes

// n:=3;
//print allinvs(OL);
//[ 1, 1, 2, 2, 1/2 ]
//print allinvs(OA);
//[ 1, 1, 1, 1, 1/24 ]

// Could only compute Left/Right classes for maximal orders
//n:=4; load "lip.m";
//print invs(OL);
//[ 1, 2, 20 ]
//print allinvs(OA);
//[ 1, 2, 2, 2, 5/48 ]

// Already here I couldn't compute the Left/Right classes for maximal orders
//n:=5; load "lip.m";
//print invs(OL);
//[ 1, 58, 2234880 ]
//print invs(OA);
//[ 1, 58, 1455/32 ]

//n:=6; load "lip.m";
//print invs(OL);
//
//print allinvs(OA);
//

