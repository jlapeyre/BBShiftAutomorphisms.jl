using Test
using LinearAlgebra: I
using BBShiftAutomorphisms
const CM = BBShiftAutomorphisms.CheckMatrices

@testset "BB codes registry" begin
    codes = list_bb_codes()
    @test codes == [(12, 6), (12, 12), (6, 6), (15, 3), (9, 6), (30, 6), (21, 18)]
    @test length(codes) == 7
end

@testset "Shift matrices (Smatrix)" begin
    s = Smatrix(5)
    eye5 = Matrix{Int}(I, 5, 5)
    @test size(s) == (5, 5)
    @test s^5 == eye5
    for k in 1:4
        @test s^k != eye5
    end
    s2 = Smatrix(6, 2)
    eye6 = Matrix{Int}(I, 6, 6)
    @test size(s2) == (6, 6)
    @test s2^3 == eye6
    @test s2 != eye6
    @test s2^2 != eye6
end

@testset "Xmatrix/Ymatrix properties" begin
    l = 3
    m = 4
    x = Xmatrix(l, m)
    y = Ymatrix(l, m)
    n = l * m
    eyen = Matrix{Int}(I, n, n)
    @test size(x) == (n, n)
    @test size(y) == (n, n)
    @test x * y == y * x
    @test x^l == eyen
    @test y^m == eyen
end

@testset "BBCodeData construction and basic checks" begin
    for (l, m) in list_bb_codes()
        data = get_bb_code_data(l, m)
        @test data isa BBCodeData
        @test get_lm(data) == (l, m)
        @test CM.check_nkd(data) === nothing
        @test CM.check_matrix_properties(data) === nothing
        (A, B) = get_AB(data)
        n = l * m
        @test size(A) == (n, n)
        @test size(B) == (n, n)
        (x, y) = get_xy(data)
        @test size(x) == (n, n)
        @test size(y) == (n, n)
    end
end

@testset "Named code access" begin
    bb1 = get_bb_code(12, 6)
    bb2 = get_bb_code(:gross)
    @test bb1 isa BBCode
    @test bb2 isa BBCode
    @test bb1.lm == (12, 6)
    @test bb2.lm == (12, 6)
    data_named = get_bb_code_data(:gross)
    @test data_named isa BBCodeData
    @test get_lm(data_named) == (12, 6)
end

@testset "Paper verification" begin
    @test verify_paper() == true
end
