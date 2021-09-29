#using Printf



RNum = Union{Integer,Rational}

"Convert x to integer if it is rational with denominator 1."
demote(x::Number) = x
demote(x::RNum) = ((denominator(x) == 1) ? numerator(x) : x)

rmul(x::Number, y::Number) = x * y
rmul(x::Integer, y::Rational) = demote(x*y)
rmul(x::Rational, y::Integer) = demote(x*y)

rdiv(x::Number, y::Number) = x / y
rdiv(x::RNum, y::RNum) = demote(x//y)



abstract type Operator end

# Default multiplication rule for operators, which can be specialised.
#
# Assumption for the moment: multiplication returns a pair
#   (c, [ops...])
# consisting of a coefficient and a possibly empty list of operators.
#
# By default we return c = 1 and a list of the same two operators given as
# inputs, i.e., we just concatenate the operators.
Base.:*(x::Operator, y::Operator) = (1, [x, y])

# This controls how lists of operators are multiplied.
# It is not very general at the moment.
function join_ops(opsx::Array{Operator,1}, opsy::Array{Operator,1})
    opx = opsx[end]
    opy = opsy[1]
    (c, opxy) = opx * opy

    if c == 0
        return (0, [])
    end

    ops = vcat(opsx[1:end-1], opxy, opsy[2:end])

    return (c, ops)
end

# Default equality test. This should be specialised for specific types of
# operators (e.g., projectors), so if this default one is called it means the
# arguments are not the same type of operator so they are not equal.
Base.:(==)(x::Operator, y::Operator) = false


# Default order. This again should be specialised for operators of the same
# type so here we just see if the type names are ordered lexicographically.
function Base.isless(x::Operator, y::Operator)
    return isless(nameof(typeof(x)), nameof(typeof(y)))
end

Base.adjoint(o::Operator) = conj(o)

Base.show(io::IO, o::Operator) = print_op(io, o)



alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z']

function party2string(p::Integer)
    base = length(alphabet)
    chars = Array{Char,1}()

    while (p > 0)
        p -= 1
        push!(chars, alphabet[1 + p % base])
        p = div(p, base)
    end

    return String(reverse!(chars))
end



struct Projector <: Operator
    output::Integer
    input::Integer
end

function print_op(io::IO, p::Projector)
    @printf io "P%d|%d" p.output p.input
end

function print_op(io::IO, p::Projector, party::Integer)
    @printf io "P%s%d|%d" party2string(party) p.output p.input
end

Base.hash(p::Projector, h::UInt) = hash((p.output, p.input), h)

function Base.:(==)(p::Projector, q::Projector)
    return (p.input == q.input) && (p.output == q.output)
end

function Base.isless(p::Projector, q::Projector)
    pi, qi = p.input, q.input
    return (pi < qi) || ((pi == qi) && (p.output < q.output))
end

function Base.:*(p::Projector, q::Projector)
    if p.input == q.input
        if p.output == q.output
            return (1, [p])
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [p, q])
    end
end

Base.conj(p::Projector) = p



struct Zbff <: Operator
    index::Integer
    conj::Bool
end

function print_op(io::IO, p::Zbff)
    @printf io "Z%s%d" (p.conj ? "*" : "") p.index
end

function print_op(io::IO, p::Zbff, party::Integer)
    @printf io "Z%s%s%d" party2string(party) (p.conj ? "*" : "") p.index
end

Base.hash(p::Operator, h::UInt) = hash((p.index, p.conj), h)

function Base.:(==)(p::Zbff, q::Zbff)
    return (p.index == q.index) && (p.conj == q.conj)
end

function Base.isless(p::Zbff, q::Zbff)
    pi, qi = p.index, q.index
    return (pi < qi) || ((pi == qi) && !p.conj && q.conj)
end

Base.conj(p::Zbff) = Zbff(p.index, !p.conj)



struct Monomial
    word::Array{Tuple{Integer,Array{Operator,1}},1}
end

function Monomial(party::Integer, operator::Operator)
    return Monomial([(party, [operator])])
end

Id = Monomial([])

Base.iterate(m::Monomial) = iterate(m.word)
Base.iterate(m::Monomial, state) = iterate(m.word, state)

Base.length(m::Monomial) = length(m.word)

Base.hash(m::Monomial, h::UInt) = hash(m.word, h)

function Base.show(io::IO, m::Monomial)
    if isempty(m)
        print(io, " Id")
    else
        for (party, ops) in m
            for o in ops
                print(io, " ")
                print_op(io, o, party)
            end
        end
    end
end

function order(m::Monomial)
    result = 0

    for (_, ops) in m.word
        result += length(ops)
    end

    return result
end



Base.:(==)(x::Number, y::Monomial) = (x == 1) && isempty(y)

Base.:(==)(x::Monomial, y::Number) = (y == 1) && isempty(x)

Base.:(==)(x::Monomial, y::Monomial) = (x.word == y.word)

function Base.isless(x::Monomial, y::Monomial)
    ox, oy = order(x), order(y)

    if ox < oy
        return true
    elseif ox == oy
        return isless(x.word, y.word)
    else
        return false
    end
end



function Base.conj(m::Monomial)
    return Monomial([(party, reverse!([conj(op) for op in ops]))
                     for (party, ops) in m])
end

function Base.adjoint(m::Monomial)
    return Monomial([(party, reverse!([adjoint(op) for op in ops]))
                     for (party, ops) in m])
end

Base.zero(m::Monomial) = Polynomial()



function operators(m::Monomial)
    return [Monomial(p, o) for (p, ops) in m for o in ops]
end


Base.:*(x::Number, y::Monomial) = Polynomial(x, y)
Base.:*(x::Monomial, y::Number) = Polynomial(y, x)

function Base.:*(x::Monomial, y::Monomial)
    coeff = 1

    if (M = length(x)) == 0
        return y
    end

    if (N = length(y)) == 0
        return x
    end

    j = 1
    k = 1

    word = Array{Tuple{Integer,Array{Operator,1}},1}()

    while (j <= M) && (k <= N)
        (px, opsx) = x.word[j]
        (py, opsy) = y.word[k]

        if px < py
            push!(word, x.word[j])
            j += 1
        elseif py < px
            push!(word, y.word[k])
            k += 1
        else
            (c, ops) = join_ops(opsx, opsy)

            if c == 0
                return 0
            end

            coeff *= c

            if !isempty(ops)
                push!(word, (px, ops))
            end

            j += 1
            k += 1
        end
    end

    append!(word, x.word[j:end])
    append!(word, y.word[k:end])

    m = Monomial(word)

    return (coeff == 1) ? m : Polynomial(m, word)
end




IndexRange = Union{UnitRange{Integer},
                   StepRange{Integer,Integer},
                   Array{Integer}}



function projector(party, output, input)
    return Monomial(party, Projector(output, input))
end

function projector(party, output::IndexRange, input::Integer)
    return [projector(party, o, input) for o in output]
end

function projector(party, output::Integer, input::IndexRange)
    return [projector(party, output, i) for i in input]
end

function projector(party, output::IndexRange, input::IndexRange)
    return [projector(party, o, i) for o in output, i in input]
end



function zbff(party, index, conj=false)
    return Monomial(party, Zbff(index, conj))
end

function zbff(party, index::IndexRange, conj=false)
    return [projector(party, a, input) for a in index]
end



struct Polynomial
    terms::Dict{Monomial,Number}
end

Polynomial() = Polynomial(Dict{Monomial,Number}())

Polynomial(x::Number) = Polynomial((x != 0) ? Dict(Id => x) : Dict())

Polynomial(x::Monomial) = Polynomial(Dict(x => 1))

function Polynomial(x::Number, y::Monomial)
    return (x != 0) ? Polynomial(Dict(y => x)) : 0
end

Polynomial(x::Polynomial) = x

Polynomial(x::Base.Generator) = Polynomial(Dict(x))

Base.iterate(x::Polynomial) = iterate(x.terms)
Base.iterate(x::Polynomial, state) = iterate(x.terms, state)

Base.getindex(x::Polynomial, y::Monomial) = get(x.terms, y, 0)

function Base.setindex!(x::Polynomial, y::Number, z::Monomial)
    if y == 0
        delete!(x.terms, z)
    else
        x.terms[z] = y
    end
end

function Base.copy(x::Polynomial)
    return Polynomial(copy(x.terms))
end

function Base.show(io::IO, p::Polynomial)
    if isempty(p)
        print(io, " 0")
    else
        for (m, c) in sort(p)
            print(io, " + (", c, ")")
            show(io, m)
        end
    end
end


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
    for (m, c) in y
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
    for (m, c) in z
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
    for (m, c) in y
        x[m] -= c
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



Base.:-(x::Monomial) = Polynomial(Dict(y => -1))
Base.:-(x::Polynomial) = Polynomial((m, -c) for (m, c) in y)



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



function Base.:*(x::Number, y::Polynomial)
    return (x != 0) ? Polynomial((m, x*c) for (m, c) in y) : 0
end

function Base.:*(x::Polynomial, y::Number)
    return (y != 0) ? Polynomial((m, y*c) for (m, c) in x) : 0
end

function Base.:*(x::Monomial, y::Polynomial)
    z = Polynomial()

    for (m, c) in y
        addmul!(z, c, x*m)
    end

    return z
end

function Base.:*(x::Polynomial, y::Monomial)
    z = Polynomial()

    for (m, c) in x
        addmul!(x, c, m*y)
    end

    return z
end

function Base.:*(x::Polynomial, y::Polynomial)
    z = Polynomial()

    for (mx, cx) in x
        for (my, cy) in y
            addmul!(z, cx*cy, mx*my)
        end
    end

    return z
end



Base.:/(x::Monomial, y::Number) = Polynomial(1/y, x)

function Base.:/(x::Polynomial, y::Number)
    divs = ((m, c/y) for (m, c) in x)
    return Polynomial((m, c) for (m, c) in divs if c != 0)
end



Base.:(==)(x::Number, y::Polynomial) = isempty(y - x)
Base.:(==)(x::Polynomial, y::Number) = isempty(x - y)

Base.:(==)(x::Monomial, y::Polynomial) = isempty(y - x)
Base.:(==)(x::Polynomial, y::Monomial) = isempty(x - y)

Base.:(==)(x::Polynomial, y::Polynomial) = isempty(x - y)



function Base.conj(x::Polynomial)
    return Polynomial((conj(m), conj(c)) for (m, c) in x)
end

function Base.adjoint(x::Polynomial)
    return Polynomial((adjoint(m), adjoint(c)) for (m, c) in x)
end

Base.zero(::Polynomial) = Polynomial()



Base.length(p::Polynomial) = length(p.terms)



Base.sort(g::Base.Generator; kws...) = sort!([x for x in g]; kws...)
Base.sort(g::Base.KeySet; kws...) = sort!([x for x in g]; kws...)

"Return pairs (m, c) of monomials and coefficients of polynomial."
function Base.sort(p::Polynomial)
    return sort!([(m, c) for (m, c) in p], by=first)
end



"Return the monomials in polynomial x."
monomials(p::Polynomial) = keys(p.terms)




function substitute!(p::Polynomial, x::Polynomial, m::Monomial)
    if ((pm = p[m]) == 0) || ((xm = x[m]) == 0)
        return p
    end

    pdivx = rdiv(pm, xm)

    for (mx, c) in x
        p[mx] -= rmul(c, pdivx)
    end

    p[m] = 0

    return p
end

substitute(p::Polynomial, x, m) = substitute!(copy(p), x, m)
