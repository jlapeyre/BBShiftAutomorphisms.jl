using Test

# Group monomorphism: Z6×Z6 → GL(6, GF(2))
@testset "verify_ax_ay_properties" begin
    Axy = get_Ax_Ay_gross()

    (ax, ay) = (Axy.Ax, Axy.Ay)
    @assert iszero(ax * ay - ay * ax)
    @assert isone(ax^6)
    @assert isone(ay^6)
    # Error will be thrown if ax, ay not invertible
    @assert isone(ax * inv(ax))
    @assert isone(ay * inv(ay))
    @assert Adelta(Axy, 1, 0) * Adelta(Axy, 0, 1) == Adelta(Axy, 1, 1)
    @assert Adelta(Axy, 0, 1) * Adelta(Axy, 1, 0) == Adelta(Axy, 1, 1)

    verify_group_operation = (Axy, a, b, c, d) ->
        Aδ(Axy, a, b) * Aδ(Axy, c, d) == Aδ(Axy, a + c, b + d)

    for (a, b, c, d) in Iterators.product(0:5, 0:5, 0:5, 0:5)
        @assert verify_group_operation(Axy, a, b, c, d)
    end

    M = Set([Aδ(Axy, a, b) for a in 0:5 for b in 0:5])
    @assert length(M) == 36

    # Kernel is trivial
    ker = [(a, b) for a in 0:5 for b in 0:5 if isone(Aδ(Axy, a, b))]
    @assert ker == [(0, 0)]
end
