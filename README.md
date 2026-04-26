# BBShiftAutomorphisms

### Automorphisms in Tour de Gross

Check the algorithm for determining how many generators are required
to implement a given shift automorphism.

This "algorithm" or formula for computing number of generators is verified as follows.

* Translate the list of six generators from Tour de gross into tuples of the exponents of
  x and y.
* Programatically add the inverses to this list.
* Generate a list of all elements of Z6xZ6 that can be represented by either one or two
  of the generators (including inverses). Make a dictionary mapping each element to the
  number of generators required.
* Verify that all elements of Z6xZ6 are represented in the result.
* Write a function that computes the number of generators required; something
  more efficient than brute force.
* Verify, for each element of Z6xZ6 that the number obtained by brute force agrees
  with the efficient function.

### Check matrices for bivariate bicycle code

Follows: High-threshold and low-overhead fault-tolerant quantum memory (Br2024)
http://arxiv.org/abs/2308.07915v2

Explicit construction of check matrices for gross and two gross codes.

Verification of some properties presented in the paper


### Testing

You can run the functions that perform verifications by hand.
They can also be run by using the standard package test.
For example
```julia
julia> using Pkg
julia> Pkg.test("BBShiftAutomorphisms")
```
