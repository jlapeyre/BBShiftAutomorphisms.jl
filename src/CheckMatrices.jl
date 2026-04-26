"""
Most, but not all, of the quantities here can be computed using type `Int` for integers,
then doing `% 2`. This does not work for `rank(hx)` and `rank(hz)`.
So we use the implementation of finite fields in Nemo.

This code follows: High-threshold and low-overhead fault-tolerant quantum memory (Bra+2024)
http://arxiv.org/abs/2308.07915v2

Data characterizing a BB code,
as well as a bit of derived data.
All data is rather small.

For gross code [[144,12,12]]  Bra+2024, Table 3
And two-gross code [[288,12,18]]
And Bra+2024 eq 1
"""
module CheckMatrices

using LinearAlgebra: I, rank
using Dictionaries: Dictionary
import Nemo

# GF2 matrix math is >10x faster with Nemo
# import AbstractAlgebra

"""
    GF2_REF :: Base.RefValue{Any}

Process-local storage for the Nemo GF(2) field object.

It is assigned in __init__() via Nemo.finite_field(2)[1] so that the
C-side pointer is valid for the current Julia session. Do not cache
the field elsewhere; access it through gf2().
"""
const GF2_REF = Base.RefValue{Any}()

"""
    gf2() -> Nemo.fmpz_mod_ctx

Return the session-valid Nemo GF(2) field stored in GF2_REF.

Nemo’s field objects wrap C pointers that change between sessions,
so the field is created in __init__() and retrieved via this
accessor to avoid stale pointers after precompilation.
"""
gf2() = GF2_REF[]

function __init__()
    return GF2_REF[] = Nemo.finite_field(2)[1]
end

# Slow
# const GF2 = AbstractAlgebra.GF(2)

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
struct BBCode{AFT, BFT}
    lm::Tuple{Int, Int} # (l, m)
    nkd::Tuple{Int, Int, Int} # [[n, k, d]]
    inv_enc_rate::Int # Inverse of Net Encoding rate
    Afunc::AFT # Eg. x^2 + y + y^2
    Bfunc::BFT
end

# struct BBCode
#     lm::Tuple{Int,Int} # (l, m)
#     nkd::Tuple{Int,Int,Int} # [[n, k, d]]
#     inv_enc_rate::Int # Inverse of Net Encoding rate
#     Afunc::Function # Eg. x^2 + y + y^2
#     Bfunc::Function
# end

"""
    BBCodeXY{XYT}

Lightweight container for the pair of shift matrices x and y associated
to a BB code. These are permutation matrices over GF(2) constructed
from the code dimensions (l, m).
"""
struct BBCodeXY{XYT}
    x::XYT
    y::XYT
end

function Base.show(io::IO, ::MIME"text/plain", xy::BBCodeXY)
    (nx, mx) = size(xy.x)
    (ny, my) = size(xy.y)
    return print(io, "BBCodeXY($nx x $mx, $ny x $my)")
end

"""
    BBCodeAB{ABT}

Lightweight container for the pair of polynomial-evaluated matrices A and B.
For a given BB code, A(x, y) and B(x, y) are formed by evaluating the
bivariate polynomials Afunc and Bfunc on the permutation matrices x and y
over GF(2).
"""
struct BBCodeAB{ABT}
    A::ABT
    B::ABT
end

_as_tuple(XY::BBCodeXY) = (XY.x, XY.y)
_as_tuple(AB::BBCodeAB) = (AB.A, AB.B)

function Base.show(io::IO, ::MIME"text/plain", ab::BBCodeAB)
    (na, ma) = size(ab.A)
    (nb, mb) = size(ab.B)
    return print(io, "BBCodeAB($na x $ma, $nb x $mb)")
end

"""
    struct BBCodeData

Store the data characterizing a BB code, together with the
associated matrices `x, y, A, B`.
"""
struct BBCodeData{ABPT, XYPT}
    bb_code::BBCode
    xy::ABPT
    AB::XYPT
end

function Base.show(io::IO, mm::MIME"text/plain", bb_mats::BBCodeData)
    (l, m) = bb_mats.bb_code.lm
    print(io, "BBCodeData((l,m)=($l, $m), ")
    show(io, mm, bb_mats.xy)
    print(io, "), ")
    show(io, mm, bb_mats.AB)
    return print(io, ")")
end

"""
    _make_x_and_y(bb_code::BBCode)

Construct the permutation matrices x and y over GF(2) for the given code.
The sizes are determined by (l, m) in bb_code.lm, with x = S_l ⊗ I_m and
y = I_l ⊗ S_m. Returns a pair (x, y).
"""
_make_x_and_y(bb_code::BBCode) = _x_and_y_matrix(bb_code.lm...)

"""
    _make_A_and_B(bb_code::BBCode, x, y)

Evaluate the code’s bivariate polynomials on the shift matrices x and y.
Returns a pair (A, B) where A = Afunc(x, y) and B = Bfunc(x, y), both
as matrices over GF(2).
"""
_make_A_and_B(bb_code::BBCode, x, y) = (bb_code.Afunc(x, y), bb_code.Bfunc(x, y))

"""
    _make_xyAB(bb_code::BBCode)

Build the tuple (x, y, A, B) for the given code. The shift matrices x and y
are formed from (l, m). The matrices A and B are then obtained by evaluating
the code polynomials at (x, y). Returns (x, y, A, B).
"""
function _make_xyAB(bb_code::BBCode)
    (x, y) = _make_x_and_y(bb_code)
    (A, B) = _make_A_and_B(bb_code, x, y)
    return (x, y, A, B)
end

BBCodeData(::Nothing) = nothing
function BBCodeData(bb_code::BBCode)
    (x, y, A, B) = _make_xyAB(bb_code)
    xy = BBCodeXY(x, y)
    AB = BBCodeAB(A, B)
    return BBCodeData{typeof(xy), typeof(AB)}(bb_code, xy, AB)
end

"""
    get_lm(code_data::BBCodeData) -> (l, m)

Return the code dimensions (l, m) from the stored BBCode.
"""
get_lm(bb_code_data::BBCodeData) = bb_code_data.bb_code.lm

"""
    get_nkd(code_data::BBCodeData) -> (n, k, d)

Return the triple [[n, k, d]] characterizing the code length, number of
logical qubits, and distance, respectively.
"""
get_nkd(bb_code_data::BBCodeData) = bb_code_data.bb_code.nkd

"""
    get_xy(code_data::BBCodeData) -> (x, y)

Return the pair of shift matrices (x, y) over GF(2) associated with the code.
These satisfy x y = y x and have orders l and m, respectively.
"""
get_xy(code_data::BBCodeData) = _as_tuple(code_data.xy)

"""
    get_AB(code_data::BBCodeData) -> (A, B)

Return the pair of code matrices (A, B) over GF(2) obtained by evaluating
the bivariate polynomials defined by the code data at (x, y).
"""
get_AB(code_data::BBCodeData) = _as_tuple(code_data.AB)

"""
    check_nkd(code_data::BBCodeData; verbose=false)

Validate the code parameters against the constructed check matrices.

- Verify `n = 2 l m`.
- Form `H^X` and `H^Z` and check `rank(H^X) = rank(H^Z)`.
- Verify `k = n - 2 rank(H^Z)`.

Throws an assertion error if any condition fails. With `verbose=true`, prints
intermediate ranks and computed `k`.
"""
function check_nkd(code_data::BBCodeData; verbose = false)
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
    return @assert k == k_calc
end

"""
    matrix_order(M)

Return smallest `n` such that `M^n == I`.

The order is computed by brute force. Assume exact arithmetic. i.e. do not use `isapprox`.
"""
function matrix_order(M)
    n_max = 100
    Mp = M
    for n in 1:n_max
        isone(Mp) && return n
        Mp = Mp * M
    end
    throw(ArgumentError(lazy"Matrix order > max $n_max"))
end

"""
    check_matrix_properties(code_data::BBCodeData; verbose=false)

Check the core algebraic relations implied by the BB construction.

Verifies:
- `x` and `y` commute, with orders `l` and `m`, respectively.
- `A` and `B` commute over GF(2), so `A B + B A = 0`.
- With `H^X = [A B]` and `H^Z = [B' A']`, confirm `H^X (H^Z)' = 0`.

Throws an assertion error if any property fails. With verbose=true, prints
selected diagnostics.
"""
function check_matrix_properties(code_data::BBCodeData; verbose = false)
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
    # In GF2, q + q is identically zero. And A and B commute. Therefore,
    @assert iszero(A * B + B * A)
    verbose && @show iszero(A * B + B * A)

    (hx, hz) = HX_HZ(code_data)
    hz_transpose = transpose(hz)
    @assert hz_transpose == vcat(B, A)

    # hx * hz' = AB + BA, and AB + BA is zero.
    return @assert iszero(hx * hz_transpose)
end

# Codes defined in (Bra+2024) Table 3, page 8
# The gross code and two-gross code are moved to the front of the list.

"""
    _bb_codes() -> Dictionary{Tuple{Int,Int}, BBCode}

Internal registry of known BB codes keyed by (l, m).

The data are taken from Bra+2024 (see Table 3 and related text). Used by get_bb_code and
list_bb_codes. The gross code and two-gross code are moved to the front of the list.
"""
function _bb_codes()
    codes = [
        # gross code
        BBCode(
            (12, 6), (144, 12, 12), 24,
            (x, y) -> x^3 + y + y^2,
            (x, y) -> y^3 + x + x^2
        ),
        # two-gross code
        BBCode(
            (12, 12), (288, 12, 18), 48,
            (x, y) -> x^3 + y^2 + y^7,
            (x, y) -> y^3 + x + x^2
        ),
        BBCode(
            (6, 6), (72, 12, 6), 12,
            (x, y) -> x^3 + y + y^2,
            (x, y) -> y^3 + x + x^2
        ),
        BBCode(
            (15, 3), (90, 8, 10), 23,
            (x, y) -> x^9 + y + y^2,
            (x, y) -> 1 + x^2 + x^7,
        ),
        BBCode(
            (9, 6), (108, 8, 10), 27,
            (x, y) -> x^3 + y + y^2,
            (x, y) -> y^3 + x + x^2
        ),
        BBCode(
            (30, 6), (360, 12, 24), 60,  # < 24, not exactly 24
            (x, y) -> x^9 + y + y^2,
            (x, y) -> y^3 + x^25 + x^26
        ),
        BBCode(
            (21, 18), (756, 16, 34), 95,  # < 34, not exactly 34
            (x, y) -> x^3 + y^10 + y^17,
            (x, y) -> y^5 + x^3 + x^19
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

"""
    get_bb_code_name(l, m) -> Symbol | nothing

Return a symbolic name for a code when one is defined. Currently returns
:gross for (12, 6) and :twogross for (12, 12). Returns nothing when no
conventional name is assigned.
"""
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
    throw(ArgumentError(lazy"Unknown BB code name '$code_name'"))
end

"""
    get_bb_code_data(l, m)::Union{BBCodeData, Nothing}
    get_bb_code_data(code_name::Symbol)::Union{BBCodeData, Nothing}

Call `get_bb_code(l, m)`, then generate the matrices A, B, x, y.
"""
get_bb_code_data(args...) = BBCodeData(get_bb_code(args...))

"""
    Smatrix(m::Integer, expt::Integer=1)::Matrix{Int}

Return an `m x m` cyclic shift matrix raised to the `expt` power.

Elements at `(i, i+1 mod m)` are equal to 1.
All others are equal to zero.

`S_l` in Bra+2024, Section 4.
"""
function Smatrix(m::Integer, expt::Integer = 1)
    s = zeros(Int, m, m)
    for i in 1:m
        j = (i + expt - 1) % m + 1
        s[i, j] = 1
    end
    return s
end

"""
    Xmatrix(l, m; expt=1) -> Matrix{Int}

Return the integer permutation matrix `x = Sₗ ⊗ Iₘ`

This is the shift on the l-index, repeated across m blocks. Intended for
sanity checks and visualization in `Z`, not over GF(2).
"""
function Xmatrix(l, m, expt = 1)
    return kron(Smatrix(l, expt), I(m))
end

"""
    Ymatrix(l, m; expt=1) -> Matrix{Int}

Return the integer permutation matrix `y = Iₗ ⊗ Sₘ`.
This is the shift on the m-index within each of the l blocks. Intended for
sanity checks and visualization in Z, not over GF(2).
"""
function Ymatrix(l, m, expt = 1)
    return kron(I(l), Smatrix(m, expt))
end

"""
    _x_and_y_matrix(l, m) -> (x, y)

Construct the pair of permutation matrices (x, y) over GF(2) with
x = S_l ⊗ I_m and y = I_l ⊗ S_m. These commute and have orders l and m,
respectively. Returned as Nemo matrices over GF(2).

`S_l` in Bra+2024, Section 4
"""
function _x_and_y_matrix(l, m)
    x = GF2mat(kron(Smatrix(l), I(m)))
    y = GF2mat(kron(I(l), Smatrix(m)))
    return (x, y)
end

"""
    GF2mat(m::Matrix)

Convert an integer matrix to a Nemo matrix over GF(2).

Entries of m are reduced modulo 2 and the size is preserved.
This is a convenience wrapper around Nemo.matrix(gf2(), m),
using the shared GF(2) instance returned by gf2().

# Examples
```jldoctest
julia> M = Smatrix(4)
julia> N = GF2mat(M)
julia> N * N == GF2mat(M * M)
true
```

See also: gf2, Smatrix, Xmatrix, Ymatrix.
"""
function GF2mat(m::Matrix)
    return Nemo.matrix(gf2(), m)
    #    return AbstractAlgebra.matrix(GF2, m)
end

"""
    HX_HZ(code_data::BBCodeData) -> (H^X, H^Z)

Form the X- and Z-type parity-check matrices from A and B as in Bra+2024 Eq. (2).

Returns:
- H^X = [A  B]
- H^Z = [B' A']
All matrices are over GF(2).
"""
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
function verify_paper(; verbose = false)
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
        check_nkd(code_data; verbose = verbose)
        check_matrix_properties(code_data; verbose = verbose)
        verbose && println()
    end
    return true
end

end # module CheckMatrices
