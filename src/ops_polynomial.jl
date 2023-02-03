struct Polynomial
    terms::Dict{Monomial,Number}
end

Polynomial() = Polynomial(Dict{Monomial,Number}())

Polynomial(x::Number) = Polynomial((x != 0) ? Dict(Id => demote(x)) : Dict())

Polynomial(x::Monomial) = Polynomial(Dict(x => 1))

function Polynomial(x::Number, y::Monomial)
    return (x != 0) ? Polynomial(Dict(y => demote(x))) : Polynomial()
end

Polynomial(x::Polynomial) = x

function Polynomial(x::Base.Generator)
    return Polynomial(Dict((m, demote(c)) for (c, m) in x))
end

new_polynomial(x) = Polynomial(x)
new_polynomial(p::Polynomial) = copy(p)



# Iteration over Polynomials.

function swap_mc(result)
    if isnothing(result)
        return nothing
    else
        ((m, c), state) = result
        return (Pair{Number,Monomial}(c, m), state)
    end
end

Base.iterate(x::Polynomial) = swap_mc(iterate(x.terms))

Base.iterate(x::Polynomial, state) = swap_mc(iterate(x.terms, state))



Base.hash(p::Polynomial, h::UInt) = hash(p.terms, h)

Base.getindex(x::Polynomial, y::Monomial) = get(x.terms, y, 0)

function Base.setindex!(x::Polynomial, y::Number, z::Monomial)
    if y == 0
        delete!(x.terms, z)
    else
        x.terms[z] = demote(y)
    end
end

function Base.copy(x::Polynomial)
    return Polynomial(copy(x.terms))
end



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



"Add y to polynomial x, modifying x."
function add!(x::Polynomial, y::Number)
    x[Id] += y
    return x
end

function add!(x::Polynomial, y::Monomial)
    x[y] += 1
    return x
end

function add!(x::Polynomial, y::Polynomial)
    for (c, m) in y
        x[m] += c
    end

    return x
end



"Add y*z to the polynomial x, modifying x. y has to be a number."
function addmul!(x::Polynomial, y::Number, z::Number)
    x[Id] += y*z
    return x
end

function addmul!(x::Polynomial, y::Number, z::Monomial)
    x[z] += y
    return x
end

function addmul!(x::Polynomial, y::Number, z::Polynomial)
    for (c, m) in z
        x[m] += c*y
    end

    return x
end



"Subtract y from polynomial x."
function sub!(x::Polynomial, y::Number)
    x[Id] -= y
    return x
end

function sub!(x::Polynomial, y::Monomial)
    x[y] -= 1
    return x
end

function sub!(x::Polynomial, y::Polynomial)
    for (c, m) in y
        x[m] -= c
    end

    return x
end



function mul!(x::Polynomial, y::Number)
    for (c, m) in x
        x[m] = c*y
    end

    return x
end



"Return a polynomial consisting of the sum of items in s."
function psum(s)
    z = Polynomial()

    for x in s
        add!(z, x)
    end

    return z
end

psum(m::Monomial) = Polynomial(m)



Base.:+(x::Number, y::Monomial) = add!(Polynomial(y), x)
Base.:+(x::Monomial, y::Number) = add!(Polynomial(x), y)

function Base.:+(x::Monomial, y::Monomial)
    return Polynomial((x != y) ? Dict(x => 1, y => 1) : Dict(x => 2))
end

Base.:+(x::Number, y::Polynomial) = add!(copy(y), x)
Base.:+(x::Polynomial, y::Number) = add!(copy(x), y)

Base.:+(x::Monomial, y::Polynomial) = add!(copy(y), x)
Base.:+(x::Polynomial, y::Monomial) = add!(copy(x), y)

Base.:+(x::Polynomial, y::Polynomial) = add!(copy(x), y)



Base.:-(x::Monomial) = Polynomial(Dict(x => -1))
Base.:-(x::Polynomial) = Polynomial((-c, m) for (c, m) in x)



Base.:-(x::Number, y::Monomial) = add!(-y, x)
Base.:-(x::Monomial, y::Number) = sub!(Polynomial(x), y)

function Base.:-(x::Monomial, y::Monomial)
    return Polynomial((x != y) ? Dict(x => 1, y => -1) : Dict())
end

Base.:-(x::Number, y::Polynomial) = add!(-y, x)
Base.:-(x::Polynomial, y::Number) = sub!(copy(x), y)

Base.:-(x::Monomial, y::Polynomial) = sub!(Polynomial(x), y)
Base.:-(x::Polynomial, y::Monomial) = sub!(copy(x), y)

Base.:-(x::Polynomial, y::Polynomial) = sub!(copy(x), y)



Base.:*(x::Number, y::Monomial) = Polynomial(x, y)
Base.:*(x::Monomial, y::Number) = Polynomial(y, x)

function Base.:*(x::Monomial, y::Monomial)
    (c, m) = join_monomials(x, y)

    return (c != 1) ? Polynomial(c, m) : m
end

function Base.:*(x::Number, y::Polynomial)
    return (x != 0) ? Polynomial((x*c, m) for (c, m) in y) : 0
end

function Base.:*(x::Polynomial, y::Number)
    return (y != 0) ? Polynomial((y*c, m) for (c, m) in x) : 0
end

function Base.:*(x::Monomial, y::Polynomial)
    z = Polynomial()

    for (c, m) in y
        addmul!(z, c, x*m)
    end

    return z
end

function Base.:*(x::Polynomial, y::Monomial)
    z = Polynomial()

    for (c, m) in x
        addmul!(z, c, m*y)
    end

    return z
end

function Base.:*(x::Polynomial, y::Polynomial)
    z = Polynomial()

    for (cx, mx) in x
        for (cy, my) in y
            addmul!(z, cx*cy, mx*my)
        end
    end

    return z
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



comm(x::Number, y::Number) = 0
comm(x::Number, y::Monomial) = 0
comm(x::Monomial, y::Number) = 0
comm(x, y) = x*y - y*x

acomm(x::Number, y::Number) = 2*rmul(x, y)
acomm(x::Number, y::Monomial) = Polynomial(2*x, y)
acomm(x::Monomial, y::Number) = Polynomail(2*y, x)
acomm(x, y) = x*y + y*x



Base.:(==)(x::Number, y::Polynomial) = isempty(y - x)
Base.:(==)(x::Polynomial, y::Number) = isempty(x - y)

Base.:(==)(x::Monomial, y::Polynomial) = isempty(y - x)
Base.:(==)(x::Polynomial, y::Monomial) = isempty(x - y)

Base.:(==)(x::Polynomial, y::Polynomial) = isempty(x - y)



function Base.conj(x::Polynomial)
    return Polynomial((conj(c), conj(m)) for (c, m) in x)
end

function Base.adjoint(x::Polynomial)
    return Polynomial((adjoint(c), adjoint(m)) for (c, m) in x)
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
