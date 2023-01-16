# Tests for Heisenberg and Weil representation.  Eventually write these as @test statements.

#Z2 = ZN(2); Z3 = ZN(3); Z4 = ZN(4); Z5 = ZN(5); Z6 = ZN(6); Z7 = ZN(7); Z8 = ZN(8);

@testset "heis" begin
    v = rand(Z8^2); w = rand(Z8^2); 

    @test heis(v)*heis(w)*heis(-v)*heis(-w) == heispairing(v,w)
    @test dagger(heis(v)) == heis(-v)
    @test dagger(heis2(v)) == heis2(-v)
end


@testset "weil" begin
    G = SL(2,Z5); g1 = rand(G); g2 = rand(G);
    w = weil_overlaps

    @test w(g1)*w(g2) == w(g1*g2)
    
    g = rand(G); j = rand(Z5^2); 
    @test w(g)*w(j)*w(g^-1) == w(g*j)
    
    @test weil_U(g)*heis(j)*weil_U(g^-1) == heis(g*j)
end

#f = j -> heisQ(j)
#matrix([(abelian_2cocycle(Z22,f)(a,b))//(abelian_2cocycle(Z22,f)(b,a)) for a in Z22, b in Z22])
