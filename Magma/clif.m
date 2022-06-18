n:=1;

G := Sp(2^n,2);
M := GModule(G);
A := AffineSplitExtension(M);

Q<i,j,k> := QuaternionAlgebra<Rationals() | -1,-1>;
L := 

//d:=8;
//G := SL(2,Integers(d));
//M := GModule(G);
//
//A := AffineSplitExtension(M);
//
//GG := SL(2,Integers(2*d));
//
//N := ChangeRing(GModule(GG),Integers(d));
//
//B, p1,p2 := AffineSplitExtension(N);

// WARNING: Magma treats module elements as row vectors
// so the ordering is reversed from the the usual semidirect product
// and the action is transposed.

//a := p2([1+d,0,0,1+d]);
//b := p2([1,d,0,1])*p1([0,d/2]);
//c := p2([1,0,d,1])*p1([d/2,0]);
//
//K := sub<B|a,b,c>;
//
//Order(a), Order(b), Order(c), #K, IsNormal(B,K);
//
//B0 := quo<B|K>;

//"Are A and B0 isomorphic?", IsIsomorphic(A,B0);
//
//"#A, #B0", #A, #B0;


//g1 := Random(GG);
//g2 := Random(GG);
//n1 := Random(N);
//n2 := Random(N);


//"Correct group action", p2(g1)*p1(n1)*p2(g2)*p1(n2) eq p2(g1*g2)*p1(n1*g2 + n2);


//J := sub<GG | [1+d,0,0,1+d], [1,d,0,1], [1,0,d,1]>;
//#J, IsNormal(GG,J), IsIsomorphic(G,quo<GG|J>);


// Gives false, so B0 is NOT isomorphic to 

//N4 := quo<M8|M8 ! [4,0],M8 ! [0,4]>;

// K:=PSL(3,4);
// H:=CyclicGroup(2);
// A:=AutomorphismGroup(K);
/* A.1 is an automorphism of order 2 */
// phi:= hom< H -> A | <H.1,A.1> >;
// G:=SemidirectProduct(K,H,phi);
