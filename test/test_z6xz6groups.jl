using Test
using BBShiftAutomorphisms

@testset "Z6: arithmetic and powers" begin
    @test Z6(-5).g == 1
    @test Z6(4) * Z6(5) == Z6(3)
    @test (Z6(2))^(-1) == inv(Z6(2))
    @test (Z6(5))^6 == Z6(0)
end

@testset "Z6: ordering and additive law for ^" begin
    @test sort([Z6(3), Z6(-4), Z6(0)]) == [Z6(0), Z6(2), Z6(3)]
    g = Z6(4)
    @test g^3 * g^(-5) == g^(-2)
end

@testset "Z6xZ6: negatives and ordering" begin
    @test Z6xZ6(-1, -7) == Z6xZ6(5, 5)
    @test Z6xZ6(4, 5) * Z6xZ6(-3, -4) == Z6xZ6(1, 1)
    @test sort([Z6xZ6(1, 0), Z6xZ6(0, 5), Z6xZ6(0, -1)]) == [Z6xZ6(0, 5), Z6xZ6(0, 5), Z6xZ6(1, 0)]
end

@testset "Z6xZ6s: basics and group laws" begin
    e = Z6xZ6(0, 0)
    @test e * e == e
    els = [Z6xZ6(a, b) for a in 0:5 for b in 0:5]
    for g in els
        @test e * g == g
        @test g * e == g
        @test g * inv(g) == e
        @test inv(g) * g == e
        @test inv(inv(g)) == g
    end
    for a in els, b in els, c in els
        @test (a * b) * c == a * (b * c)
    end
end

@testset "Z6xZ6s: constructors and ordering" begin
    @test Z6(7).g == 1
    @test Z6xZ6((1, 7)) == Z6xZ6(1, 1)
    expected = [Z6xZ6(a, b) for a in 0:5 for b in 0:5]
    @test sort(copy(expected)) == expected
end

@testset "Z6xZ6s: generating set properties" begin
    gens = generating_set
    @test length(gens) == 12
    pairs = [(g.x.g, g.y.g) for g in gens]
    @test length(unique(pairs)) == 12
    for g in gens
        @test inv(g) in gens
    end
end

@testset "Z6xZ6s: enumeration via generators" begin
    dict = compute_nr_generators()
    @test length(dict) == 36
    all_els = [Z6xZ6(a, b) for a in 0:5 for b in 0:5]
    @test sort(collect(keys(dict))) == all_els
    @test all(v in (1, 2) for v in values(dict))
    @test dict[Z6xZ6(0, 0)] == 2
    @test all_elements_present()
    @test all_elements_present(collect(keys(dict)))
end

@testset "Z6xZ6s: lookup for number of generators" begin
    @test check_nr_generators_lookup()
    @test nr_generators_lookup(Z6xZ6(0, 0)) == 2
    for xy in ((1, 0), (0, 1), (5, 0), (0, 5))
        @test nr_generators_lookup(Z6xZ6(xy)) == 1
    end
    for xy in ((3, 3), (0, 3), (3, 0))
        @test nr_generators_lookup(Z6xZ6(xy)) == 2
    end
    for a in 0:5
        xy = (3, a)
        if !(xy in ((3, 3), (3, 0)))
            @test nr_generators_lookup(Z6xZ6(xy)) == 1
        end
        xy = (a, 3)
        if !(xy in ((3, 3), (0, 3)))
            @test nr_generators_lookup(Z6xZ6(xy)) == 1
        end
    end
end
