PartiedOpList = Vector{Tuple{Integer,Vector{Operator}}}

struct Monomial
    word::PartiedOpList
end

function Monomial(party::Integer, operator::Operator)
    @assert party > 0
    return Monomial([(party, [operator])])
end

Monomial(party, operator::Operator) = Monomial(party_num(party), operator)

if !@isdefined(Id)
    const Id = Monomial([])
end

isidentity(m::Monomial) = isempty(m)

Base.iterate(m::Monomial) = iterate(m.word)
Base.iterate(m::Monomial, state) = iterate(m.word, state)

Base.length(m::Monomial) = length(m.word)

Base.hash(m::Monomial, h::UInt) = hash(m.word, h)

function show_nonempty_pol(io::IO, word::PartiedOpList)
    sep = ""

    for (party, ops) in word
        for o in ops
            print(io, sep)
            print(io, string(o, party))
            sep = " "
        end
    end
end

function show_pol(io::IO, word::PartiedOpList)
    if isempty(word)
        print(io, "Id")
    else
        show_nonempty_pol(io, word)
    end
end

Base.show(io::IO, m::Monomial) = show_pol(io, m.word)

degree(x::Number) = !iszero(x) ? 0 : -Inf

function degree(m::Monomial)
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
    ox, oy = degree(x), degree(y)

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


function conj_pol(x::PartiedOpList)
    return [(party, reverse!([conj(op) for op in ops]))
            for (party, ops) in x]
end

Base.conj(m::Monomial) = Monomial(conj_pol(m.word))

function Base.adjoint(m::Monomial)
    return Monomial([(party, reverse!([adjoint(op) for op in ops]))
                     for (party, ops) in m])
end

Base.zero(m::Monomial) = 0



conj_min(x::Number) = real(x)

function conj_min(m::Monomial)
    return min(m, conj(m))
end



"""
Concatenate two monomials. This is used later to decide what the result
of multiplying two monomials is.
"""
function join_pol(x::PartiedOpList, y::PartiedOpList)
    coeff = 1

    if (M = length(x)) == 0
        return (coeff, y)
    end

    if (N = length(y)) == 0
        return (coeff, x)
    end

    j = 1
    k = 1

    word = PartiedOpList()

    while (j <= M) && (k <= N)
        (px, opsx) = x[j]
        (py, opsy) = y[k]

        if px < py
            push!(word, x[j])
            j += 1
        elseif py < px
            push!(word, y[k])
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

    append!(word, x[j:end])
    append!(word, y[k:end])

    return (coeff, word)
end
