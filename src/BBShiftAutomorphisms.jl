module BBShiftAutomorphisms

using Reexport

include("Z6xZ6Groups.jl")
@reexport using .Z6xZ6Groups: Z6Group, Z6xZ6Group, compute_nr_generators, all_elements_present, generating_set,
    nr_generators_lookup, check_nr_generators_lookup

include("CheckMatrixes.jl")
@reexport using .CheckMatrixes: verify_paper, BBCode, BBCodeData, get_bb_code, list_bb_codes, get_xy, get_AB,
    get_lm,
    get_bb_code_data,
    Smatrix,
    Xmatrix,
    Ymatrix

end # module BBShiftAutomorphisms
