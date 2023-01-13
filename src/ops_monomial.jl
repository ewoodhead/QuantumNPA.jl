OpVector = Vector{Tuple{PartyVec,Vector{Operator}}}

struct Monomial
    word::OpVector
end

function Monomial(party, operator::Operator)
    return Monomial([(party_vec(party), [operator])])
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
                print(io, string(o, party))
                sep = " "
            end
        end
    end
end

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



Base.zero(m::Monomial) = 0



conj_min(x::Number) = real(x)

conj_min(m::Monomial) = min(m, conj(m))



"""
Test if vectors u and v are the same up to index n. Returns 0 if they are,
otherwise returns the first index where they differ.
"""
function same_upto(u, v, n)
    k = 0
    
    while ((k += 1) <= n)
        if u[k] != v[k]
            return k
        end
    end

    return 0
end

"""
Test if two party vecs have any parties in common, optionally starting from
given indices j0 and k0.
"""
function parties_isect(u::PartyVec,
                       v::PartyVec,
                       j=1,
                       k=1,
                       m=length(u),
                       n=length(v))
    if (j > m) || (k > n)
        return true
    end

    p = u[j]
    q = v[k]

    while true
        if p == q
            return true
        elseif p < q
            if (j += 1) > m
                return false
            end
         
            p = u[j]
        else
            if (k += 1) > n
                return false
            end

            q = v[k]
        end
    end
end

"""
Test if two party vecs should be swapped (return :swap), merged
(return :join), or left in the order they already are (return :leave).
"""
function swap_or_join(u::PartyVec, v::PartyVec)
    if (m = length(u)) == 0
        return isempty(v) ? :join : :leave
    elseif (n = length(v)) == 0
        return :swap
    end

    if u[1] < v[1]
        return :leave
    end

    if (m == n) && ((k0 = same_upto(u, v, m)) == 0)
        return :join
    else
        k0 = 1
    end

    return parties_isect(u, v, k0, k0, m, n) ? :leave : :swap
end

function insert_at(x::OpVector, j0::Int, m::Int, q::PartyVec)
    for j in m:-1:j0
        p = x[j][1]
        result = swap_or_join(p, q)

        if result === :join
            return (j, :join)
        elseif result === :leave
            return (j+1, :leave)
        end
    end

    return (j0, :leave)
end

function join_words(x::OpVector, y::OpVector)
    coeff = 1

    if (m = length(x)) == 0
        return (coeff, y)
    elseif (n = length(y)) == 0
        return (coeff, x)
    end

    word = copy(x)
    j = 1

    for (k, qv) in enumerate(y)
        (q, v) = qv
        (j, action) = insert_at(word, j, m, q)

        if action === :join
            (c, w) = join_ops(word[j][2], v)

            if iszero(c)
                return (0, OpVector[])
            end

            coeff *= c

            if !isempty(w)
                word[j] = (q, w)
                j += 1
            else
                deleteat!(word, j)
                j = 1
                m -= 1
            end
        else
            insert!(word, j, qv)
            j += 1
            m += 1
        end
    end

    return (coeff, word)
end

"""
Concatenate two monomials. This is used later to decide what the result
of multiplying two monomials is.
"""
function join_monomials(x::Monomial, y::Monomial)
    (c, word) = join_words(x.word, y.word)
    return (c, Monomial(word))
end



function word_adjoint(w::OpVector)
    
end

Base.adjoint(m::Monomial) = Monomial(word_adjoint(m.word))

Base.conj(m::Monomial) = Base.adjoint(m::Monomial)



function ctrace(m::Monomial)
    coeff = 1

    pcops = [(p, trace(ops)) for (p, ops) in m.word]

    for (_, (c, _)) in pcops
        coeff = rmul(coeff, c)
    end

    m = Monomial([(p, ops) for (p, (_, ops)) in pcops])

    return (coeff == 1) ? m : (coeff, m)
end
