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

# This needs to be redefined for operators that aren't Hermitian.
Base.conj(o::Operator) = o

Base.adjoint(o::Operator) = conj(o)

Base.show(io::IO, o::Operator) = print(io, string(o))



alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z']



"Return string representation of a party, e.g., `party_str(3) = \"C\"`."
function party_str(p::Integer)
    base = length(alphabet)
    chars = Array{Char,1}()

    while (p > 0)
        p -= 1
        push!(chars, alphabet[1 + p % base])
        p = div(p, base)
    end

    return String(reverse!(chars))
end

party_str(s::String) = s
party_str(s::Symbol) = string(s)


pos(c::Char) = findfirst(isequal(c), alphabet)

"Return integer representation of a party, e.g., `party_num(\"C\") = 3`."
function party_num(s::String)
    base = length(alphabet)
    party = 0

    for c in s
        party = pos(c) + base * party
    end

    return party
end

party_num(c::Char) = pos(c)
party_num(s::Symbol) = party_num(string(s))
party_num(n::Integer) = n

"Split string into party and rest, e.g. \"AB1\" -> (\"AB\", \"1\")."
function split_party(s::String)
    k = findlast(in(alphabet), s)
    return (s[1:k], s[k+1:end])
end

split_party(s::Symbol) = split_party(string(s))

function range_expr(expr)
    
end



struct Monomial
    word::Array{Tuple{Integer,Array{Operator,1}},1}
end

function Monomial(party::Integer, operator::Operator)
    @assert party > 0
    return Monomial([(party, [operator])])
end

Monomial(party, operator::Operator) = Monomial(party_num(party), operator)

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
                print(io, string(o, party))
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

Base.zero(m::Monomial) = 0



IndexRange = Union{UnitRange{<:Integer},Array{<:Integer}}

getfields(expr) = (expr.head == :tuple) ? expr.args : [expr]
getfieldnames(fields) = map(argname, fields)

argname(arg::Symbol) = arg
argname(arg::Expr) = arg.args[1]

function conjfalse(field)
    if argname(field) === :conj
        return Expr(:kw, field, false)
    else
        return field
    end
end

instance_fields(instance, names) = [:($instance.$f) for f in names]

function fmt_remove(fmt, s::Symbol)
    return Expr(fmt.head, filter(!isequal(s), fmt.args)...)
end

function parse_fmt(fmt)
    if fmt.head === :string
        return (fmt, fmt_remove(fmt, :party))
    else
        fmts = fmt.args
        return (fmts[1], fmts[2])
    end
end

function stringfdef(name, fmt)
    fieldnames = [f for f in fmt.args if f isa Symbol]

    args = if (:party in fieldnames)
               (:(x::$name), :(party::Integer))
           else
               (:(x::$name),)
           end

    function fix_special(f)
        if f === :conj
            return :((x.$f ? "*" : ""))
        elseif f === :party
            return :(party_str(party))
        else
            return :(x.$f)
        end
    end
    
    bindings = (:($f = $(fix_special(f))) for f in fieldnames)

    return :(function Base.string($(args...))
                 $(bindings...)
                 return $fmt
             end)
end

"""
Define a new type of operator with a given name (e.g., Projector), fields,
and format string. In addition to generating the struct defintion this also
generates method definitions for the following generic functions:

  * Base.hash,
  * Base.:(==),
  * Base.isless,
  * Base.string,

as well as a constructor method with the name in lowercase (e.g., projector)
that creates a monomial containing a single operator associated to a given
party.

If one of the fields is named conj a Base.conj method is also generated.
"""
macro operator(ctor::Expr, fmt, order=nothing)
    name = ctor.args[1]
    fields = ctor.args[2:end]
    lcname = Symbol(lowercase(string(name)))

    fieldnames = getfieldnames(fields)

    if :conj in fieldnames
        cxfields = replace(instance_fields(:x, fieldnames),
                           :(x.conj) => :(!x.conj))
        conjf = :( Base.conj(x::$name) = $name($(cxfields...)) )
        cargs = map(conjfalse, fields)
    else
        conjf = nothing
        cargs = fields
    end

    order = (isnothing(order) ? reverse(fieldnames) : getfields(order))
    xfields = :($(instance_fields(:x, order)...),)
    yfields = :($(instance_fields(:y, order)...),)

    (fmt_party, fmt_noparty) = parse_fmt(fmt)
    strf_party = stringfdef(name, fmt_party)
    strf_noparty = stringfdef(name, fmt_noparty)

    methods = filter(!isnothing, (strf_party, strf_noparty, conjf))

    return quote
        struct $name <: Operator
            $(fields...)
        end
        Base.hash(x::$name, h::UInt) = hash($xfields, h)
        Base.:(==)(x::$name, y::$name) = ($xfields == $yfields)
        Base.isless(x::$name, y::$name) = ($xfields < $yfields)
        $(methods...)
        function $(esc(lcname))(party, $(cargs...))
            return Monomial(party, $name($(fieldnames...)))
        end
    end
end



@operator Dichotomic(input::Integer) ("$party$input", "/$input")

function Base.:*(x::Dichotomic, y::Dichotomic)
    return (x.input == y.input) ? (1, []) : (1, [x, y])
end

dichotomic(party, input::IndexRange) = [dichotomic(party, z) for z in input]

function parse_dichotomic(expr)
    if expr isa Symbol
        (party, input) = split_party(expr)
        name = Symbol(party)
        return :($(esc(expr)) = dichotomic($party, $(parse(Int, input))))
    else
        name = expr.args[1]
        range = expr.args[2]
        return :($(esc(name)) = dichotomic($(QuoteNode(name)), $range))
    end
end

macro dichotomic(expr...)
    exprs = map(parse_dichotomic, expr)
    return quote $(exprs...) end
end



@operator(Fourier(input::Integer, power::Integer, d::Integer),
          "$party$input^$power",
          (d, input, power))

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



@operator(Projector(output::Integer, input::Integer),
          "P$party$output|$input")

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

function projector(party, output::IndexRange, input::Integer;
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

function projector(party, output::IndexRange, input::IndexRange;
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



@operator(KetBra(outputl::Integer, outputr::Integer, input::Integer),
          "|$outputl><$outputr|$input$party")

function Base.:*(p::Projector, k::KetBra)
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

function Base.:*(k::KetBra, p::Projector)
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

function Base.:*(k::KetBra, l::KetBra)
    input = k.input

    if input == l.input
        kor, lol = k.outputr, l.outputl
        
        if kor == lol
            kol, lor = k.outputl, l.outputr

            if kol != lor
                return (1, [KetBra(kol, lor, input)])
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

Base.conj(k::KetBra) = KetBra(k.outputr, k.outputl, k.input)

function ketbra(party, outputl::Integer, outputr::Integer, input::Integer)
    op = ((outputl != outputr) ?
          KetBra(outputl, outputr, input) :
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



@operator Unitary(index::Integer, conj::Bool) "U$party$conj$index"

function Base.:*(u::Unitary, v::Unitary)
    if (u.index != v.index) || (u.conj == v.conj)
        return (1, [u, v])
    else
        return (1, [])
    end
end

function unitary(party, index::IndexRange, conj=false)
    return [unitary(party, i, conj) for i in index]
end



@operator Zbff(index::Integer, conj::Bool) "Z$party$conj$index"

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
    return (x != 0) ? Polynomial(Dict(y => demote(x))) : Polynomial()
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

    return (coeff == 1) ? m : Polynomial(coeff, m)
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



function add_monomials!(table::Dict{Integer,Set{Monomial}},
                        itr)
    for x in itr
        add_monomials!(table, x)
    end

    return table
end

function add_monomials!(table::Dict{Integer,Set{Monomial}},
                        p::Polynomial)
    for m in monomials(p)
        add_monomials!(table, m)
    end

    return table
end

"Update table of parties -> monomials with operators in a given monomial."
function add_monomials!(table::Dict{Integer,Set{Monomial}},
                        m::Monomial)
    for (p, ops) in m
        if !haskey(table, p)
            table[p] = Set{Monomial}()
        end

        for o in ops
            push!(table[p], Monomial(p, o))
        end
    end

    return table
end



"Return all the individual (order 1) operators in the arguments."
operators(s...; by_party::Bool=false) = operators(s; by_party=by_party)

function operators(itr; by_party::Bool=false)
    if !by_party
        return Set{Monomial}(flatten(map(operators, itr)))
    else
        return add_monomials!(Dict{Integer,Set{Monomial}}(), itr)
    end
end

function operators(p::Polynomial; by_party::Bool=false)
    if !by_party
        return Set{Monomial}(flatten(operators(m) for m in monomials(p)))
    else
        return add_monomials!(Dict{Integer,Set{Monomial}}(), p)
    end
end

"Return all the individual operators making up a monomial."
function operators(m::Monomial; by_party::Bool=false)
    if !by_party
        return Set{Monomial}(Monomial(p, o)
                             for (p, ops) in m for o in ops)
    else
        return add_monomials!(Dict{Integer,Set{Monomial}}(), m)
    end
end

function operators(x::Number; by_party::Bool=false)
    return !by_party ? Set{Monomial}() : Dict{Integer,Set{Monomial}}()
end

function operators(table::Dict{Integer,Set{Monomial}}; by_party::Bool=false)
    if !by_party
        return reduce(union!, values(table); init=Set{Monomial}())
    else
        return table
    end
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
    PA = projector(1, 1:d, 1:2; full=true)
    PB = projector(2, 1:d, 1:2; full=true)

    d1 = d - 1

    function p(x, y, k)
        return psum(PA[1+a, x] * PB[1+mod(a+k,d), y] for a in 0:d1)
    end

    return psum(mul!(p(1,1,k) + p(2,1,-k-1) + p(2,2,k) + p(1,2,-k)
                     - p(1,1,-k-1) - p(2,1,k) - p(2,2,-k-1) - p(1,2,k+1),
                     (1 - 2*k//d1))
                for k in 0:(div(d,2)-1))
end
