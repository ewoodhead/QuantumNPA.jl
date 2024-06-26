module QuantumNPA

export npa_min, npa_max, npa2sdp, sdp2jump, npa2jump,
    npa_moment, npa_level,
    direct_sum,
    set_solver!,
    Id, @dichotomic,
    dichotomic, fourier, unitary, projector, zbff,
    Monomial, Polynomial, monomials, coefficients, operators,
    cglmp

using Base.Iterators #  flatten, zip (filter?)
using BlockDiagonals
using Combinatorics #  powerset
using JuMP
using LinearAlgebra
using SCS
using SparseArrays

include("operators.jl")
include("ops_predefined.jl")
include("npa.jl")

end
