module LogicalShifts

using ..CheckMatrices: GF2mat

# Equation (37) in Tour de gross. Page 47
#
const A_x_gross = [
    0 1 0 1 0 0;
    0 1 0 0 0 1;
    0 0 1 1 0 0;
    1 1 0 1 1 0;
    0 1 0 0 1 0;
    1 1 1 1 0 1
]

const A_y_gross = [
    1 0 0 0 0 1;
    1 1 1 0 0 1;
    0 0 0 0 1 0;
    0 1 0 0 0 0;
    0 1 1 0 0 1;
    0 0 1 1 0 1;
]

struct AxAy{T}
    Ax::T
    Ay::T
end

get_Ax_Ay_gross() = AxAy(GF2mat(A_x_gross), GF2mat(A_y_gross))

Adelta(Axy, a, b) = Axy.Ax^mod(a, 6) * Axy.Ay^mod(b, 6)
const Aδ = Adelta

# Group monomorphism: Z6×Z6 → GL(6, GF(2))
function verify_ax_ay_properties(Axy)
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

    ker = [(a, b) for a in 0:5 for b in 0:5 if isone(Aδ(Axy, a, b))]
    return @assert ker == [(0, 0)]
end

end # module LogicalShifts
