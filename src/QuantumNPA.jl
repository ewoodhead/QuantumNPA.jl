module QuantumNPA

export @dichotomic, Id, Monomial, Polynomial, cglmp, coefficients, comm,
    dichotomic, direct_sum, fourier, isidentity, monomials, npa2jump,
    npa2sdp, npa_level, npa_max, npa_min, npa_moment, operators, projector,
    sdp2jump, set_solver!, unitary, generic

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
