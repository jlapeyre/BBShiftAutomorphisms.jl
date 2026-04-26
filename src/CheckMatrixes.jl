module CheckMatrixes

using LinearAlgebra: I, rank
using Dictionaries: Dictionary
using Nemo: finite_field, matrix, identity_matrix

# Most, but not all, of the quantities here can be computed using type `Int` for intgers,
# then doing `% 2`. This does not work for `rank(hx)` and `rank(hz)`.
# So we use the implementation of finite fields in Nemo

# This code follows: High-threshold and low-overhead fault-tolerant quantum memory (Bra+2024)
# http://arxiv.org/abs/2308.07915v2

# Data characterizing a BB code,
# as well as a bit of derived data.
# All data is rather small.
#
# For gross code [[144,12,12]]  Bra+2024, Table 3
# And two gross code [[288,12,18]]
# And Bra+2024 eq 1

"""
    struct BBCode

Data characterizing a BB code, as well as a bit of derived data.

Nothing stored here is resource intensive (neither CPU, nor memory).
Associated matrices may be computed and stored in `struct BBCodeData`.
Data included here is copied from Table 3 in (Bra+24).

## Fields
- lm::Tuple{Int,Int} : `(l, m)`
- nkd::Tuple{Int,Int,Int} : `[[n, k, d]]`
- inv_enc_rate::Int : Inverse of Net Encoding rate
- Afunc, Bfunc : functions mapping (x, y) to matrices A and B. These are polynomials in x and y.
"""
struct BBCode
    lm::Tuple{Int,Int} # (l, m)
    nkd::Tuple{Int,Int,Int} # [[n, k, d]]
    inv_enc_rate::Int # Inverse of Net Encoding rate
    Afunc::Function # Eg. x^2 + y + y^2
    Bfunc::Function
end

_make_x_and_y(bb_code::BBCode) = _x_and_y_matrix(bb_code.lm...)
_make_A_and_B(bb_code::BBCode, x, y) = (bb_code.Afunc(x, y), bb_code.Bfunc(x, y))

function _make_xyAB(bb_code::BBCode)
    (x, y) = _make_x_and_y(bb_code)
    (A, B) = _make_A_and_B(bb_code, x, y)
    return (x, y, A, B)
end

struct BBCodeXY{XYT}
    x::XYT
    y::XYT
end

function Base.show(io::IO, ::MIME"text/plain", xy::BBCodeXY)
    (nx, mx) = size(xy.x)
    (ny, my) = size(xy.y)
    print(io, "BBCodeXY($nx x $mx, $ny x $my)")
end

struct BBCodeAB{ABT}
    A::ABT
    B::ABT
end

_as_tuple(XY::BBCodeXY) = (XY.x, XY.y)
_as_tuple(AB::BBCodeAB) = (AB.A, AB.B)

function Base.show(io::IO, ::MIME"text/plain", ab::BBCodeAB)
    (na, ma) = size(ab.A)
    (nb, mb) = size(ab.B)
    print(io, "BBCodeAB($na x $ma, $nb x $mb)")
end

"""
    struct BBCodeData

Store the data characterizing a BB code, together with the
associated matrices `x, y, A, B`.
"""
struct BBCodeData{ABT, XYT}
    bb_code::BBCode
    xy::BBCodeXY
    AB::BBCodeAB
end

function Base.show(io::IO, mm::MIME"text/plain", bb_mats::BBCodeData)
    (l, m) = bb_mats.bb_code.lm
    print(io, "BBCodeData((l,m)=($l, $m), ")
    show(io, mm, bb_mats.xy)
    print(io, "), ")
    show(io, mm, bb_mats.AB)
    print(io, ")")
end

BBCodeData(::Nothing) = nothing
function BBCodeData(bb_code::BBCode)
    (x, y, A, B) = _make_xyAB(bb_code)
    xy = BBCodeXY(x, y)
    AB = BBCodeAB(A, B)
    return BBCodeData{typeof(AB), typeof(xy)}(bb_code, xy, AB)
end

get_lm(bb_code_data::BBCodeData) = bb_code_data.bb_code.lm
get_nkd(bb_code_data::BBCodeData) = bb_code_data.bb_code.nkd
get_xy(code_data::BBCodeData) = _as_tuple(code_data.xy)
get_AB(code_data::BBCodeData) = _as_tuple(code_data.AB)

function check_nkd(code_data::BBCodeData; verbose=false)
    (l, m) = get_lm(code_data)
    (n, k, d) = get_nkd(code_data)
    n_calc = 2 * l * m
    @assert n == n_calc
    (hx, hz) = HX_HZ(code_data)
    rank_hz = rank(hz)
    rank_hx = rank(hx)
    verbose && @show (rank_hz, rank_hx)
    @assert rank_hx == rank_hz
    k_calc = n - 2 * rank_hz
    verbose && @show (k_calc = n - 2 * rank_hz, k)
    @assert k == k_calc
end

"""
    matrix_order(M)

Return smallest `n` such that `M^n == I`.

Assume exact arithemetic. i.e. do not use `isapprox`.
"""
function matrix_order(M)
    n_max = 100
    Mp = M
    for n in 1:n_max
        isone(Mp) && return n
        Mp = Mp * M
    end
    error(lazy"Matrix order > max $n_max")
end

function check_matrix_properties(code_data::BBCodeData; verbose=false)
    (x, y) = get_xy(code_data)

    # x = Sₗ ⊗ Iₘ, y = Iₗ ⊗ Sₘ => x and y commute
    @assert x * y == y * x

    (l, m) = get_lm(code_data)
    @assert matrix_order(x) == l
    @assert matrix_order(y) == m
    verbose && @show (matrix_order(x), matrix_order(y))

    (A, B) = get_AB(code_data)
    # A and B commute
    @assert A * B == B * A
    # In GF2, q + q is identicaly zero. And A and B commute. Therefore,
    @assert iszero(A * B + B * A)
    verbose && @show iszero(A * B + B * A)

    (hx, hz) = HX_HZ(code_data)
    hz_transpose = transpose(hz)
    @assert hz_transpose == vcat(B, A)

    # hx * hz' = AB + BA, and AB + BA is zero.
    @assert iszero(hx * hz_transpose)
end

# Codes defined in (Bra+2024) Table 3, page 8
# The gross code and twogross code are moved to the front of the list.
function _bb_codes()
    codes = [
        # gross code
        BBCode((12, 6), (144, 12, 12), 24,
               (x,y) -> x^3 + y + y^2,
               (x,y) -> y^3 + x + x^2
               ),
        # twogross code
        BBCode((12, 12), (288, 12, 18), 48,
               (x,y) -> x^3 + y^2 + y^7,
               (x,y) -> y^3 + x + x^2
               ),
        BBCode((6, 6), (72, 12, 6), 12,
               (x,y) -> x^3 + y + y^2,
               (x,y) -> y^3 + x + x^2
               ),
        BBCode((15, 3), (90, 8, 10), 23,
               (x,y) -> x^9 + y + y^2,
               (x,y) -> 1 + x^2 + x^7,
               ),
        BBCode((9, 6), (108, 8, 10), 27,
               (x,y) -> x^3 + y + y^2,
               (x,y) -> y^3 + x + x^2
               ),
        BBCode((30, 6), (360, 12, 24), 60,  # < 24, not exactly 24
               (x,y) -> x^9 + y + y^2,
               (x,y) -> y^3 + x^25 + x^26
               ),
        BBCode((21, 18), (756, 16, 34), 95,  # < 34, not exactly 34
               (x,y) -> x^3 + y^10 + y^17,
               (x,y) -> y^5 + x^3 + x^19
               ),
    ]
    return Dictionary([code.lm for code in codes], codes)
end

const _BB_CODES = _bb_codes()

"""
    list_bb_codes()

List `(l, m)` tuples for BB codes, for which we have data.
"""
list_bb_codes() = collect(keys(_BB_CODES))

"""
    get_bb_code(l, m)::Union{BBCode, Nothing}

Return data characterizing BB code for `(l, m)`.

This relies on the fact that `(l, m)` uniquely identifies a code
in our data set.
"""
get_bb_code(l, m) = get(_BB_CODES, (l, m), nothing)

function get_bb_code_name(l, m)
    (l, m) == (12, 6) && return :gross
    (l, m) == (12, 12) && return :twogross
    return nothing
end

"""
    get_bb_code(code_name::Symbol)

Get code by name: either `:gross` or `:twogross`.
"""
function get_bb_code(code_name::Symbol)
    code_name === :gross && return get_bb_code(12, 6)
    code_name === :twogross && return get_bb_code(12, 12)
    error(lazy"Unknown BB code name '$code_name'")
end

"""
    get_bb_code_data(l, m)::Union{BBCodeData, Nothing}
    get_bb_code_data(code_name::Symbol)::Union{BBCodeData, Nothing}

Call `get_bb_code(l, m)`, then generate the matrices A, B, x, y.
"""
get_bb_code_data(args...) = BBCodeData(get_bb_code(args...))

# S_l in Bra+2024, Section 4
# Return Matrix{Int}, not matrix over GF2

"""
    Smatrix(m::Integer, expt::Integer=1)::Matrix{Int}

Return an `m x m` cyclic shift matrix raised to the `expt` power
Elements at `(i, i+1 mod m)` are equal to 1.
All others are equal to zero.
"""
function Smatrix(m::Integer, expt::Integer=1)
    s = zeros(Int, m, m)
    for i in 1:m
        j = (i + expt - 1) % m + 1
        s[i, j] = 1
    end
    return s
end

# x = Sₗ ⊗ Iₘ
function Xmatrix(l, m, expt=1)
    kron(Smatrix(l, expt), I(m))
end

# y = Iₗ ⊗ Sₘ
function Ymatrix(l, m, expt=1)
    kron(I(l), Smatrix(m, expt))
end

# S_l in Bra+2024, Section 4
# Return matrices over GF2
function _x_and_y_matrix(l, m)
    GF2 = finite_field(2)[1]
    x = matrix(GF2, kron(Smatrix(l), I(m)))
    y = matrix(GF2, kron(I(l), Smatrix(m)))
    return (x, y)
end

# Check matrices
# Bra+2024 eq 2, H^X, H^Z
function HX_HZ(bb_code_data::BBCodeData)
    (A, B) = get_AB(bb_code_data)
    hx = hcat(A, B)
    # M' is adjoint in Julia. But this is not defined for Nemo: finite_field
    # Must use `transpose`.
    hz = hcat(transpose(B), transpose(A))
    return (hx, hz)
end

"""
    verify_paper(; verbose=false)

Verify that I understand what the math in Bra+2024 means.

Verify certain assertions made in the paper. Do this for each
of several BB codes, including most in Table 3.

See Sec. 4 Bivariate Bicycle quantum LDPC codes
"""
function verify_paper(; verbose=false)
    for (l, m) in list_bb_codes()
        bb_code = get_bb_code(l, m)
        if verbose
            bb_code_name = get_bb_code_name(l, m)
            !isnothing(bb_code_name) && @info "BB Code name: $bb_code_name"
            @info "(l, m) = ($l, $m)"
            (n, k, d) = bb_code.nkd
            @info "[[n, k, d]] = [[$n, $k, $d]]"
        end
        code_data = BBCodeData(bb_code)
        check_nkd(code_data; verbose=verbose)
        check_matrix_properties(code_data; verbose=verbose)
        verbose && println()
    end
    return true
end

end # module CheckMatrixes
