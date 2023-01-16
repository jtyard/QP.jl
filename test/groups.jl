using QP 
using Oscar


@testset "groups" begin

    @test is_isomorphic(SL(2,2),symmetric_group(3))

    @test is_isomorphic(PSL(2,3),alternating_group(4))

    @test is_isomorphic(PSL(2,4),symmetric_group(4))

    @test is_isomorphic(PSL(2,5),alternating_group(5))

    @test is_isomorphic(PGL(2,5),symmetric_group(5))

end

#[dim(V) for V in irreducible_modules(SL(2,ZN(4)))]