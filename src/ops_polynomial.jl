Scalar = Number
Coefficient = Union{Scalar,AbstractVector,AbstractMatrix}



struct Polynomial
    cfsize::Union{Tuple{},Tuple{Int},Tuple{Int,Int}}
    blockstruct::Vector{Tuple{Int,Int}}
    terms::Dict{Monomial,Coefficient}
end



# Determine blocksizes from cfsize

blockstruct(::Tuple{}) = []
blockstruct((x,)::Tuple{Int}) = [(x, x)]
blockstruct(x::Tuple{Int,Int}) = [x]



# Determine cfsize from blockstruct

cfsize(blockstruct::Vector{Tuple{Int,Int}}) = reduce(.+, blockstruct)



# Alternative constructors.

# Constructors with subsets of the struct fields (cfsize, blockstruct,
# terms).

function Polynomial(cfsize::Tuple, terms::Dict{Monomial,Coefficient})
    return Polynomial(cfsize, blockstruct(cfsize), terms)
end

function Polynomial(blockstruct::Vector{Tuple{Int,Int}},
                    terms::Dict{Monomial,Coefficient})
    return Polynomial(cfsize(blockstruct), blockstruct, terms)
end

function Polynomial(cfsize::Tuple, blockstruct::Vector{Tuple{Int,Int}})
    return Polynomial(cfsize,
                      blockstruct,
                      Dict{Coefficient,Monomial}())
end

Polynomial(cfsize::Tuple=()) = Polynomial(cfsize, blockstruct(cfsize))

function Polynomial(blockstruct::Vector{Tuple{Int,Int}})
    return Polynomial(cfsize(blockstruct), blockstruct)
end



# Constructors to create polynomials out of coefficients and monomials.

function Polynomial(x::Coefficient)
    cfsize = size(x)
    terms = (!iszero(x) ? Dict(Id => demote(x)) : Dict())
    return Polynomial(cfsize, blockstruct(cfsize), terms)
end

function Polynomial(x::BlockDiagonal)
    terms = (!iszero(x) ? Dict(Id => demote(x)) : Dict())
    return Polynomial(size(x), blockstruct(x), terms)
end

function Polynomial(m::Monomial)
    return Polynomial((), Dict(m => 1))
end

function Polynomial(c::Coefficient, m::Monomial)
    if !iszero(c)
        return Polynomial(size(c),
                          Dict{Monomial,Coefficient}(m => demote(c)))
    else
        return Polynomial(size(c))
    end
end

Polynomial(p::Polynomial) = p

function Polynomial(x::Base.Generator; cfsize::Tuple=())
    result = Iterators.peel(x)

    if isnothing(result)
        return Polynomial(cfsize)
    else
        ((c, m), rest) = result
        sz = size(c)
        terms = Dict(m => demote(c))

        for (c, m) in rest
            terms[m] = demote(c)
        end

        return Polynomial(sz, terms)
    end
end

function Polynomial(x::Base.Generator, p::Polynomial)
    return Polynomial(x, cfsize=p.cfsize)
end



# Accessors.

Base.size(p::Polynomial) = p.cfsize
blocksizes(p::Polynomial) = p.blockstruct



# Construct a new polynomial out of arg that can be safely modified.

new_polynomial(x) = Polynomial(x)
new_polynomial(p::Polynomial) = copy(p)



# Iteration over Polynomials.

function swap_mc(result)
    if isnothing(result)
        return nothing
    else
        ((m, c), state) = result
        return ((c, m), state)
    end
end

Base.iterate(x::Polynomial) = swap_mc(iterate(x.terms))
Base.iterate(x::Polynomial, state) = swap_mc(iterate(x.terms, state))



Base.hash(p::Polynomial, h::UInt) = hash(p.terms, h)



Base.iszero(p::Polynomial) = isempty(p.terms)


# Make it possible to access/assign polynomial coefficients via [] accessor.

function zero_coeff(p::Polynomial)
    cfsize = p.cfsize
    
    if cfsize === ()
        return 0
    else
        return zeros(Int, cfsize)
    end
end

function Base.getindex(p::Polynomial, m)
    terms = p.terms

    if haskey(terms, m)
        return terms[m]
    else
        return zero_coeff(p)
    end
end

function Base.setindex!(p::Polynomial, c, m)
    if iszero(c)
        delete!(p.terms, m)
    else
        p.terms[m] = demote(c)
    end
end



# Copy a polynomial. We need to copy the dictionary of terms since these are
# often mutated, but not the others.                      
                      
function Base.copy(p::Polynomial)
    return Polynomial(p.cfsize, p.blockstruct, copy(p.terms))
end



# Default printing of polynomials at the REPL.

function Base.show(io::IO, p::Polynomial)
    if isempty(p)
        print(io, "0")
    else
        c2s = firstcoeff2string

        for (c, m) in sort(p)
            print(io, c2s(c))
            show(io, m)
            c2s = coeff2string
        end
    end
end



"""
Degree of a polynomial.

degree(0) returns negative infinity. With this definition the following rules
hold even if one or both of P or Q are zero:

  degree(P + Q) == max(degree(P), degree(Q))
  degree(P - Q) <= max(degree(P), degree(Q))
  degree(P * Q) == degree(P) + degree(Q)
"""
function degree(p::Polynomial)
    return !iszero(p) ? maximum(degree.(monomials(p))) : -Inf
end

degree_less(n::Integer) = (x -> (degree(x) < n))



"Delete entry m in terms if c is zero, otherwise set terms[m] = c."
function set_nonzero!(terms::Dict, c, m)
    if iszero(c)
        delete!(terms, m)
    else
        terms[m] = c
    end
end
    


# Low level functions to add to a polynomial.

function add!(p::Polynomial, c::Coefficient, m::Monomial)
    @assert size(c) === p.cfsize
    return add_unsafe!(p, c, m)
end

function add_unsafe!(p::Polynomial, c::Coefficient, m::Monomial)
    if iszero(c)
        return p
    end

    terms = p.terms

    if haskey(terms, m)
        c += terms[m]
    end

    set_nonzero!(terms, c, m)

    return p
end

function add!(p::Polynomial, c::Coefficient, q::Polynomial)
    @assert p.cfsize === q.cfsize

    for (cq, m) in q
        add_unsafe!(p, cq*c, m)
    end

    return p
end

function add!(p::Polynomial, q::Polynomial)
    @assert p.cfsize === q.cfsize

    for (c, m) in q
        add_unsafe!(p, c, m)
    end

    return p
end



# Low level functions to subtract from a polynomial.

function sub!(p::Polynomial, c, m)
    @assert size(c) === p.cfsize
    return sub_unsafe!(p, c, m)
end

function sub_unsafe!(p::Polynomial, c, m)
    if iszero(c)
        return p
    end

    terms = p.terms

    if haskey(terms, m)
        c = terms[m] - c
    else
        c = -c
    end

    set_nonzero!(terms, c, m)

    return p
end

function sub!(p::Polynomial, x, q::Polynomial)
    @assert p.cfsize === q.cfsize

    for (c, m) in q
        sub_unsafe!(p, c*x, m)
    end

    return p
end

function sub!(p::Polynomial, q::Polynomial)
    @assert p.cfsize === q.cfsize

    for (c, m) in q
        sub_unsafe!(p, c, m)
    end

    return p
end



# Multiply a polynomial by a coefficient, modifying the polynomial.
# Multiplication of coefficients may not be commutative so we need separate
# left and right multiplication functions.

function mul!(p::Polynomial, c::Coefficient)
    if !iszero(c)
        terms = p.terms

        for (m, cp) in terms
            set_nonzero!(terms, cp*c, m)
        end
    else
        empty!(p.terms)
    end

    return p
end

function mul!(c::Coefficient, p::Polynomial)
    if !iszero(c)
        terms = p.terms

        for (m, cp) in terms
            set_nonzero!(terms, c*cp, m)
        end
    else
        empty!(p.terms)
    end

    return p
end



"Return a polynomial consisting of the sum of items in s."
function psum(s)
    if isempty(s)
        return Polynomial()
    end

    (x, rest) = Iterators.peel(s)
    z = new_polynomial(x)

    for x in rest
        add!(z, x)
    end

    return z
end

psum(m::Monomial) = Polynomial(m)



# Binary addition.

function Base.:+(x::Monomial, y::Monomial)
    return Polynomial((), (x != y) ? Dict(x => 1, y => 1) : Dict(x => 2))
end

Base.:+(m::Monomial, p::Polynomial) = add!(copy(p), 1, m)
Base.:+(p::Polynomial, m::Monomial) = add!(copy(p), 1, m)

Base.:+(x::Polynomial, y::Polynomial) = add!(copy(x), y)



# Unary negation.

Base.:-(m::Monomial) = Polynomial(-1, m)
Base.:-(p::Polynomial) = Polynomial(((-c, m) for (c, m) in p), p)



# Binary subtraction.

function Base.:-(x::Monomial, y::Monomial)
    terms = ((x != y) ? Dict(x => 1, y => -1) : Dict())
    return Polynomial((), terms)
end

Base.:-(m::Monomial, p::Polynomial) = sub!(Polynomial(m), p)
Base.:-(p::Polynomial, m::Monomial) = sub!(copy(p), 1, m)

Base.:-(p::Polynomial, q::Polynomial) = sub!(copy(p), q)



# Multiplication.

Base.:*(c::Coefficient, m::Monomial) = Polynomial(c, m)
Base.:*(m::Monomial, c::Coefficient) = Polynomial(c, m)

function Base.:*(x::Monomial, y::Monomial)
    (c, m) = join_monomials(x, y)
    return (c != 1) ? Polynomial(c, m) : m
end

function Base.:*(c::Coefficient, p::Polynomial)
    if iszero(c)
        return 0
    else
        pairs = ((c*cp, m) for (cp, m) in p)
        return Polynomial((c, m) for (c, m) in pairs if !iszero(c))
    end
end

function Base.:*(p::Polynomial, c::Coefficient)
    if iszero(c)
        return 0
    else
        pairs = ((cp*c, m) for (cp, m) in p)
        return Polynomial((c, m) for (c, m) in pairs if !iszero(c))
    end
end

function Base.:*(m::Monomial, p::Polynomial)
    q = Polynomial(p.cfsize)

    for (c, mp) in p
        add!(q, c, m*mp)
    end

    return q
end

function Base.:*(p::Polynomial, m::Monomial)
    q = Polynomial(p.cfsize)

    for (c, mp) in p
        add!(q, c, mp*m)
    end

    return q
end

cf_mul_size(::Tuple{}, ::Tuple{}) = ()
cf_mul_size(::Tuple{}, x) = x
cf_mul_size(x, ::Tuple{}) = x

function cf_mul_size((m,)::Tuple{Int}, (l, n)::Tuple{Int, Int})
    @assert l == 1
    return (m, n)
end

function cf_mul_size((m, l)::Tuple{Int,Int}, (n,)::Tuple{Int})
    @assert l == 1
    return (m, n)
end

function cf_mul_size((m, k)::Tuple{Int,Int}, (l, n)::Tuple{Int,Int})
    @assert k == l
    return (m, n)
end



function Base.:*(p::Polynomial, q::Polynomial)
    r = Polynomial(cf_mul_size(p.cfsize, q.cfsize))

    for (cp, mp) in p
        for (cq, mq) in q
            add!(r, cp*cq, mp*mq)
        end
    end

    return r
end



Base.:/(x::Monomial, y::Number) = Polynomial(rdiv(1, y), x)

function Base.:/(x::Polynomial, y::Number)
    divs = ((rdiv(c, y), m) for (c, m) in x)
    return Polynomial((c, m) for (c, m) in divs if c != 0)
end



function Base.:^(x::Union{Monomial,Polynomial}, p::Integer)
    @assert p >= 0

    result = ((p % 2 == 1) ? x : Id)
    p >>= 1

    while p > 0
        x *= x

        if p % 2 == 1
            result *= x
        end

        p >>= 1
    end

    return result
end



comm(::Scalar, ::Scalar) = 0
comm(::Scalar, ::Monomial) = 0
comm(::Monomial, ::Scalar) = 0
comm(x, y) = x*y - y*x

acomm(x::Scalar, y::Scalar) = 2*rmul(x, y)
acomm(x::Scalar, m::Monomial) = Polynomial(2*x, m)
acomm(m::Monomial, x::Scalar) = Polynomial(2*x, m)
acomm(x, y) = x*y + y*x



Base.:(==)(m::Monomial, p::Polynomial) = iszero(p - m)
Base.:(==)(p::Polynomial, m::Monomial) = iszero(p - m)

Base.:(==)(p::Polynomial, q::Polynomial) = iszero(p - q)



function Base.conj(p::Polynomial)
    return Polynomial((conj(c), conj(m)) for (c, m) in p)
end

function Base.adjoint(p::Polynomial)
    return Polynomial((adjoint(c), adjoint(m)) for (c, m) in p)
end

Base.zero(::Polynomial) = Polynomial()

Base.length(p::Polynomial) = length(p.terms)



function conj_min(p::Polynomial)
    return psum(conj_min(c) * conj_min(m) for (c, m) in p)
end



function trace(m::Monomial)
    (c, tm) = ctrace(m)
    return (c == 1) ? tm : Polynomial(c, tm)
end

trace(p::Polynomial) = psum(c*trace(m) for (c, m) in p)
