using BBShiftAutomorphisms
using Test

@testset "BBShiftAutomorphisms.jl" begin
    # Check matrices for BB code
    @test verify_paper()

    # Test by enumeration that all elements of Z6xZ6 are implemented by no more
    # than two generators.
    @test all_elements_present()

    # Check that the lookup function for number of generators agrees with
    # what was computed by brute force.
    @test check_nr_generators_lookup()
end
