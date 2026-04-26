# BBShiftAutomorphisms

BBShiftAutomorphisms is a small Julia package used to validate the algebra and closed‑form rules behind bivariate bicycle (BB) codes. I wrote it to cross‑check and debug a Rust implementation that counts the minimum number of generators needed to realize a given sequence of operations.

## Installation

```julia
using Pkg
Pkg.develop(url = "https://github.com/jlapeyre/BBShiftAutomorphisms.jl")
```

### Automorphisms in Tour de gross

Check the algorithm for determining how many generators are required
to implement a given shift automorphism.

This "algorithm" (a closed‑form rule for the number of generators) is verified as follows.

* Translate the six generators from Tour de gross into `(x, y)` exponent pairs.
* Add their inverses programmatically to the list of generators.
* Enumerate all elements of Z₆×Z₆ reachable by one or two generators; map each element to the minimal count.
* Verify all elements of Z₆×Z₆ are covered.
* Implement conditional logic that returns the number of generators.
* Cross‑check the logic against the brute‑force enumeration.

## Check matrices for bivariate bicycle codes (Bra+24)

Verify assertions made in  (Bra+24).

* Build shift matrices `x = S_l ⊗ I_m` and `y = I_l ⊗ S_m` over GF(2).
* Evaluate the code polynomials `A(x, y)` and `B(x, y)` (from Bra+24, Table 3).
* Form `H^X = [A  B]` and `H^Z = [B'  A']` and verify:
  * `x y = y x`; orders are `l` and `m`.
  * `A B = B A` and, over GF(2), `A B + B A = 0`.
  * `H^X (H^Z)' = 0`.
  * `n = 2 l m` and `k = n − 2 rank(H^Z)`.
* Verify some invariants of the group homomorphism (monomorphism) from Z₆×Z₆ → GL(6, GF(2)) in Eq. (37) in TdG

## Quickstart

```julia
using BBShiftAutomorphisms

list_bb_codes()
verify_paper(; verbose = true)
get_bb_code_data(:gross) |> get_AB
```

### Testing

You can run the functions that perform verifications by hand.
They can also be run by using the standard package test.
For example
```julia
julia> using Pkg
julia> Pkg.test("BBShiftAutomorphisms")
```

### Notes:

* GF(2) arithmetic: we use Nemo for matrices over GF(2) so linear‑algebra routines (e.g. `rank`) run over the correct field without manual `% 2`.
* AbstractAlgebra is lighter and pure Julia; Nemo (via FLINT C) makes tests ~20× faster here, so we use Nemo.

### References
* Tour de gross http://arxiv.org/abs/2506.03094v1
* High-threshold and low-overhead fault-tolerant quantum memory (Bra+2024) http://arxiv.org/abs/2308.07915v2
