# Tests for Heisenberg and Weil representation.  Eventually write these as @test statements.

Z2 = ZN(2); Z3 = ZN(3); Z4 = ZN(4); Z5 = ZN(5); Z6 = ZN(6); Z7 = ZN(7); Z8 = ZN(8);


N = 4;

N2 = Int(N//2);

V = Z8^2;
v = rand(V); w = rand(V); 

heis(v)*heis(w)*heis(-v)*heis(-w) == heispairing(v,w), dagger(heis(v)) == heis(-v), dagger(heis2(v)) == heis2(-v), v, det(heis(v)), det(heis2(v))

