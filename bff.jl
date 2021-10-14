#using Printf
#  @printf
#using Base.Iterators
#  flatten, zip



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
# Assumption: inputs opsx and opsy both contain at least one element.
function join_ops(opsx::Array{Operator,1}, opsy::Array{Operator,1})
    j = length(opsx)
    k = 1
    K = 1 + length(opsy)
    c = 1

    while true
        opx, opy = opsx[j], opsy[k]
        (c1, op) = opx * opy

        if c1 == 0
            return (0, [])
        end

        c *= c1
        j -= 1
        k += 1

        if (op != []) || (j == 0) || (k == K)
            ops = vcat(opsx[1:j], op, opsy[k:end])
            return (c, ops)
        end
    end
end

# Default equality test. This should be specialised for specific types of
# operators (e.g., projectors), so if this default one is called it means the
# arguments are not the same type and are therefore not equal.
Base.:(==)(x::Operator, y::Operator) = false


# Default order. This again should be specialised for operators of the same
# type so here we just see if the type names are ordered lexicographically.
function Base.isless(x::Operator, y::Operator)
    return isless(nameof(typeof(x)), nameof(typeof(y)))
end

Base.adjoint(o::Operator) = conj(o)

Base.show(io::IO, o::Operator) = print_op(io, o)



abstract type HermitianOperator <: Operator end

Base.conj(h::HermitianOperator) = h



alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z']



"Converty party number to string, e.g., `party2string(3) = 'C'`."
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

pos(c::Char) = first(indexin(c, alphabet))

"Convert string to party number, e.g., `string2party(\"C\") = 3`."
function string2party(s::String)
    base = length(alphabet)
    party = 0

    for c in s
        party = pos(c) + base * party
    end

    return party
end



struct Dichotomic <: HermitianOperator
    input::Integer
end

function print_op(io::IO, x::Dichotomic)
    @printf io "_%d" x.input
end

function print_op(io::IO, x::Dichotomic, party::Integer)
    @printf io "%s%d" party2string(party) x.input
end

Base.hash(x::Dichotomic, h::UInt) = hash(x.input, h)

Base.:(==)(x::Dichotomic, y::Dichotomic) = (x.input == y.input)

Base.isless(x::Dichotomic, y::Dichotomic) = (x.input < y.input)

function Base.:*(x::Dichotomic, y::Dichotomic)
    return (x.input == y.input) ? (1, []) : (1, [x, y])
end



struct Fourier <: Operator
    input::Integer
    power::Integer
    d::Integer
end

function print_op(io::IO, x::Fourier)
    @printf io "d%d^%d" x.input x.power
end

function print_op(io::IO, x::Fourier, party::Integer)
    @printf io "%s%d^%d" party2string(party) x.input x.power
end

Base.hash(x::Fourier, h::UInt) = hash((x.input, x.power, x.d), h)

function Base.:(==)(x::Fourier, y::Fourier)
    return (x.d == y.d) && (x.input == y.input) && (x.power == y.power)
end

function Base.isless(x::Fourier, y::Fourier)
    if x.d != y.d
        return x.d < y.d
    else
        xi = x.input
        yi = y.input
        return (xi < yi) || ((xi == yi) && (x.power < y.power))
    end
end

function Base.:*(x::Fourier, y::Fourier)
    input = x.input
    d = x.d

    if (y.d != d) || (y.input != input)
        return (1, [x, y])
    else
        p = (x.power + y.power) % d
        return (1, ((p != 0) ? [Fourier(input, p, d)] : []))
    end
end

Base.conj(x::Fourier) = Fourier(x.input, x.d - x.power, x.d)



struct Projector <: HermitianOperator
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



struct Ketbra <: Operator
    outputl::Integer
    outputr::Integer
    input::Integer
end

function print_op(io::IO, k::Ketbra)
    @printf io "|%d><%d|%d" k.outputl k.outputr k.input
end

function print_op(io::IO, k::Ketbra, party::Integer)
    @printf io "|%d><%d|%d%s" k.outputl k.outputr k.input party2string(party)
end

Base.hash(k::Ketbra, h::UInt) = hash((k.outputl, k.outputr, k.input), h)

function Base.:(==)(k::Ketbra, l::Ketbra)
    return ((k.input == l.input) && (k.outputl == l.outputr)
            && (k.outputr == l.outputr))
end

function Base.isless(k::Ketbra, l::Ketbra)
    ki, li = k.input, l.input

    return (k.input, k.outputl, k.outputr) < (l.input, l.outputl. l.outputr)
end

function Base.:*(p::Projector, k::Ketbra)
    if p.input == k.input
        if p.output == k.outputl
            return (1, [k])
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [p, k])
    end
end

function Base.:*(k::Ketbra, p::Projector)
    if k.input == p.input
        if k.outputr == p.output
            return (1, [k])
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [k, p])
    end
end

function Base.:*(k::Ketbra, l::Ketbra)
    input = k.input

    if input == l.input
        kor, lol = k.outputr, l.outputl
        
        if kor == lol
            kol, lor = k.outputl, l.outputr

            if kol != lor
                return (1, [Ketbra(kol, lor, input)])
            else
                return (1, [Projector(kol, input)])
            end
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [k, p])
    end
end

Base.conj(k::Ketbra) = Ketbra(k.outputr, k.outputl, k.input)



struct Unitary <: Operator
    index::Integer
    conj::Bool
end

function print_op(io::IO, u::Unitary)
    @printf io "U%s%d" (u.conj ? "*" : "") u.index
end

function print_op(io::IO, u::Unitary, party::Integer)
    @printf io "U%s%s%d" party2string(party) (u.conj ? "*" : "") u.index
end

Base.hash(u::Unitary, h::UInt) = hash((u.index, u.conj), h)

function Base.:(==)(u::Unitary, v::Unitary)
    return (u.index == v.index) && (u.conj == v.conj)
end

function Base.isless(u::Unitary, v::Unitary)
    ui, vi = u.index, v.index
    return (ui < vi) || ((ui == vi) && !u.conj && v.conj)
end

function Base.:*(u::Unitary, v::Unitary)
    if (u.index != v.index) || (u.conj == v.conj)
        return (1, [u, v])
    else
        return (1, [])
    end
end

Base.conj(u::Unitary) = Unitary(u.index, !u.conj)



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
    @assert party > 0
    return Monomial([(party, [operator])])
end

Id = Monomial([])

isidentity(m::Monomial) = isempty(m)

Base.iterate(m::Monomial) = iterate(m.word)
Base.iterate(m::Monomial, state) = iterate(m.word, state)

Base.length(m::Monomial) = length(m.word)

Base.hash(m::Monomial, h::UInt) = hash(m.word, h)

function Base.show(io::IO, m::Monomial)
    if isidentity(m)
        print(io, "Id")
    else
        sep = ""

        for (party, ops) in m
            for o in ops
                print(io, sep)
                print_op(io, o, party)
                sep = " "
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

    if ox != oy
        return ox < oy
    end

    for ((p1, ops1), (p2, ops2)) in zip(x, y)
        if p1 != p2
            return p1 < p2
        end

        l1, l2 = length(ops1), length(ops2)

        if l1 != l2
            return l1 > l2
        end

        for (o1, o2) in zip(ops1, ops2)
            if o1 != o2
                return o1 < o2
            end
        end
    end

    return false
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



IndexRange = Union{UnitRange{<:Integer},
                   StepRange{<:Integer,<:Integer},
                   Array{<:Integer}}



function dichotomic(party, input::Integer)
    return Monomial(party, Dichotomic(input))
end

function dichotomic(party, input::IndexRange)
    return [dichotomic(party, z) for z in input]
end



function fourier(party, input::Integer, power::Integer, d::Integer)
    @assert d > 0
    p = power % d
    return (p != 0) ? Monomial(party, Fourier(input, p, d)) : Id
end

function fourier(party, input::IndexRange, power::Integer, d::Integer)
    return [fourier(party, z, power, d) for z in input]
end

function fourier(party, input::Integer, power::IndexRange, d::Integer)
    return [fourier(party, input, power, d) for p in power]
end

function fourier(party, input::IndexRange, power::IndexRange, d::Integer)
    return [fourier(party, z, p, d) for z in input, p in power]
end



function projector(party, output::Integer, input::Integer)
    return Monomial(party, Projector(output, input))
end

function projector(party, output::IndexRange, input::Integer,
                   full::Bool=false)
    if full
        outputs = output[1:end-1]
        ps = [projector(party, o, input) for o in outputs]
        pl = Polynomial(Id)

        for p in ps
            pl[p] = -1
        end

        return vcat(ps, [pl])
    else
        return [projector(party, o, input) for o in output]
    end
end

function projector(party, output::Integer, input::IndexRange)
    return [projector(party, output, i) for i in input]
end

function projector(party, output::IndexRange, input::IndexRange,
                   full::Bool=false)
    if full
        outputs = output[1:end-1]
        ps = [projector(party, o, i) for o in outputs, i in input]
        pl = reshape([Polynomial(Id) for _ in input], 1, length(input))

        for i in input
            for p in ps[:,i]
                pl[i][p] = -1
            end
        end

        return vcat(ps, pl)
    else
        return [projector(party, o, i) for o in output, i in input]
    end
end



function ketbra(party, outputl::Integer, outputr::Integer, input::Integer)
    op = ((outputl != outputr) ?
          Ketbra(outputl, outputr, input) :
          Projector(outputl, input))
    return Monomial(party, op)
end

function ketbra(party, outputl::IndexRange, outputr::IndexRange,
                input::Integer)
    return [ketbra(party, ol, or, input) for ol in outputl, or in outputr]
end

function ketbra(party, outputl::Integer, outputr::Integer, input::IndexRange)
    return [ketbra(party, outputl, outputr, i) for i in input]
end

function ketbra(party, outputl::IndexRange, outputr::IndexRange,
                input::IndexRange)
    return [ketbra(party, ol, or, i)
            for ol in outputl, or in outputr, i in input]
end



function unitary(party, index::Integer, conj=false)
    return Monomial(party, Unitary(index, conj))
end

function unitary(party, index::IndexRange, conj=false)
    return [unitary(party, i, conj) for i in index]
end



function zbff(party, index::Integer, conj=false)
    return Monomial(party, Zbff(index, conj))
end

function zbff(party, index::IndexRange, conj=false)
    return [zbff(party, i, conj) for i in index]
end



struct Polynomial
    terms::Dict{Monomial,Number}
end

Polynomial() = Polynomial(Dict{Monomial,Number}())

Polynomial(x::Number) = Polynomial((x != 0) ? Dict(Id => demote(x)) : Dict())

Polynomial(x::Monomial) = Polynomial(Dict(x => 1))

function Polynomial(x::Number, y::Monomial)
    return (x != 0) ? Polynomial(Dict(y => demote(x))) : 0
end

Polynomial(x::Polynomial) = x

function Polynomial(x::Base.Generator)
    return Polynomial(Dict((m, demote(c)) for (m, c) in x))
end

Base.iterate(x::Polynomial) = iterate(x.terms)
Base.iterate(x::Polynomial, state) = iterate(x.terms, state)

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


num2str(x::Real) = "$x"

function num2str(x::Rational)
    a, b = numerator(x), denominator(x)

    return (b != 1) ? "$a/$b" : "$a"
end

function csgn(x::Real, p::String = "+", m::String = "-")
    return (x >= 0) ? p : m
end

function sgnnum(x::Number, p::String = "+", m::String = "-")
    xr = real(x)
    xi = imag(x)
    
    if xi == 0
        return (csgn(xr, p, m), num2str(abs(xr)))
    elseif xr == 0
        return (csgn(xi, p, m), "$(num2str(abs(xi)))im")
    else
        xis = num2str(abs(xi))

        if xr >= 0
            xrs = num2str(xr)
            s = csgn(xi)
            
            return (p, "($xrs $s $(xis)im)")
        else
            xrs = num2str(-xr)
            s = csgn(-xi)

            return (m, "($xrs $s $(xis)im)")
        end
    end
end

function firstcoeff2string(x::Number)
    if x == 1
        return ""
    elseif x == -1
        return "-"
    else
        (s, xs) = sgnnum(x, "", "-")
        return "$s$xs "
    end
end

function coeff2string(x::Number)
    if x == 1
        return " + "
    elseif x == -1
        return " - "
    else
        (s, xs) = sgnnum(x)
        return " $s $xs "
    end
end

function Base.show(io::IO, p::Polynomial)
    if isempty(p)
        print(io, "0")
    else
        c2s = firstcoeff2string

        for (m, c) in sort(p)
            print(io, c2s(c))
            show(io, m)
            c2s = coeff2string
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



function mul!(x::Polynomial, y::Number)
    for (m, c) in x
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



Base.:-(x::Monomial) = Polynomial(Dict(y => -1))
Base.:-(x::Polynomial) = Polynomial((m, -c) for (m, c) in x)



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
        addmul!(z, c, m*y)
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



Base.:/(x::Monomial, y::Number) = Polynomial(rdiv(1, y), x)

function Base.:/(x::Polynomial, y::Number)
    divs = ((m, rdiv(c, y)) for (m, c) in x)
    return Polynomial((m, c) for (m, c) in divs if c != 0)
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



Sortable = Union{Base.Generator,
                 Base.Set,
                 Base.KeySet,
                 Base.Iterators.Flatten}

Base.sort(g::Sortable; kws...) = sort!([x for x in g]; kws...)

"Return pairs (m, c) of monomials and coefficients of polynomial in order."
function Base.sort(p::Polynomial)
    return sort!([(m, c) for (m, c) in p], by=first)
end



"Return all monomials in arguments including duplicates."
all_monomials(s...) = all_monomials(s)
all_monomials(itr) = flatten(map(all_monomials, itr))
all_monomials(p::Polynomial) = keys(p.terms)
all_monomials(m::Monomial) = (m,)



"Return all the monomials in the arguments."
monomials(s...) = monomials(s)

"Return all the monomials in iterable itr."
monomials(itr) = Set(flatten(map(monomials, itr)))

"Return the monomials in polynomial x."
monomials(p::Polynomial) = keys(p.terms)

monomials(m::Monomial) = (m,)



"Return all the individual (order 1) operators in the arguments."
operators(s...) = operators(s)

operators(itr) = Set(flatten(map(operators, itr)))

function operators(p::Polynomial)
    return Set(flatten(operators(m) for m in monomials(p)))
end

"Return all the individual operators making up a monomial."
function operators(m::Monomial)
    return (Monomial(p, o) for (p, ops) in m for o in ops)
end



"Eliminate monomial m from p assuming <x> = 0. Modifies p."
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

"Remove lexicographically highest monomial in x from p assuming x = 0"
function substitute!(p::Polynomial, x::Polynomial)
    
end

substitute(p::Polynomial, x, m) = substitute!(copy(p), x, m)



function cglmp(d::Integer)
    PA = projector(1, 1:d, 1:2, true)
    PB = projector(2, 1:d, 1:2, true)

    d1 = d - 1

    function p(x, y, k)
        return psum(PA[1+a, x] * PB[1+mod(a+k,d), y] for a in 0:d1)
    end

    return psum(mul!(p(1,1,k) + p(2,1,-k-1) + p(2,2,k) + p(1,2,-k)
                     - p(1,1,-k-1) - p(2,1,k) - p(2,2,-k-1) - p(1,2,k+1),
                     (1 - 2*k//d1))
                for k in 0:(div(d,2)-1))
end
