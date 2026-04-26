module BBShiftAutomorphisms

using Reexport

include("Z6Groups.jl")
@reexport using .Z6Groups: Z6, Z6xZ6, compute_nr_generators, all_elements_present, generating_set,
    nr_generators_lookup, check_nr_generators_lookup

include("CheckMatrices.jl")
@reexport using .CheckMatrices: verify_paper, BBCode, BBCodeData, get_bb_code, list_bb_codes, get_xy, get_AB,
    get_lm,
    get_bb_code_data,
    Smatrix,
    Xmatrix,
    Ymatrix,
    gf2,
    GF2mat

include("LogicalShifts.jl")
@reexport using .LogicalShifts: A_x_gross, A_y_gross, get_Ax_Ay_gross, verify_ax_ay_properties,
    AxAy, Adelta, Aδ

end # module BBShiftAutomorphisms
