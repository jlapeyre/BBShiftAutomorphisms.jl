using BBShiftAutomorphisms
using Aqua: Aqua
import Nemo

const PACKAGE_NAME = BBShiftAutomorphisms

@testset "aqua test ambiguities" begin
    @test isempty(Test.detect_ambiguities(PACKAGE_NAME))
#    Aqua.test_ambiguities([PACKAGE_NAME, Core, Base, Nemo])
end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(PACKAGE_NAME)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(PACKAGE_NAME)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(PACKAGE_NAME)
end

@testset "aqua stale deps" begin
    Aqua.test_stale_deps(PACKAGE_NAME)
end

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(PACKAGE_NAME)
end

@testset "aqua piracies" begin
    Aqua.test_piracies(PACKAGE_NAME)
end

# @testset "aqua undocumented_names" begin
#     Aqua.test_undocumented_names(PACKAGE_NAME)
# end

# Quite slow. Unlikely we violate this in any case
# @testset "aqua persistent tasks" begin
#     Aqua.test_persistent_tasks(PACKAGE_NAME)
# end

