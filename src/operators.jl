# Multiplication and division such that x*y or x/y gives a rational if
# x and y are, or an integer if the denominator of the result is 1.
RNum = Union{Integer,Rational}

"Convert x to integer if it is rational with denominator 1."
demote(x::Number) = x
demote(x::Rational) = ((denominator(x) == 1) ? numerator(x) : x)

rmul(x::Number, y::Number) = x * y
rmul(x::Integer, y::Rational) = demote(x*y)
rmul(x::Rational, y::Integer) = demote(x*y)

rdiv(x::Number, y::Number) = x / y
rdiv(x::RNum, y::RNum) = demote(x//y)



# Definition of abstract Operator type and infrastructure for primitive
# operators.
include("ops_primitive.jl")

# Convert party numbers to strings and vice versa.
include("ops_misc.jl")

# Definition of monomials (products of operators divided into parties) and
# supporting functions.
include("ops_monomial.jl")

# Definition of scalars other than number: trace of operators.
include("ops_scalar.jl")

# Definition of polynomials (linear combinations of monomials multiplied by
# coefficients) and supporting functions.
include("ops_polynomial.jl")

# Utility functions for dealing with monomials and polynomials, e.g.,
# return all the individual monomials in a monomial or polynomial,
# reduce a list of polynomial constraints to a canonical form.
include("ops_utilities.jl")
