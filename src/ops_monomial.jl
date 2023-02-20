OpVector = Vector{Tuple{PartyVec,Vector{Operator}}}

struct Monomial
    word::OpVector
end

Id_word = OpVector()
Id = Monomial(Id_word)

function Monomial(party, operator::Operator)
    return Monomial([(party_vec(party), [operator])])
end

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


"""
Return -1, 0, or 1 depending on whether p predeces, equals, or follows q
lexicographically, assuming they are the same length.
"""
function cmp_samelen(p::Vector, q::Vector)
    for (x, y) in zip(p, q)
        if (x != y)
            return (x < y) ? -1 : 1
        end
    end

    return 0
end



"""
Return -1, 0, or 1 depending on whether p precedes, equals, or follows q.

\"Precede\" here means either:
  - p is shorter than q, or
  - p is the same length as p, and p lexicographically precedes q.
"""
function cmp_vec(p::Vector, q::Vector)
    m = length(p)
    n = length(q)

    if (m != n)
        return (m < n) ? -1 : 1
    end

    return cmp_samelen(p, q)
end

"""
Takes two pairs p_ops = (p, [x...]) and q_ops = (q, [y...]) of parties and
operators associated to those parties and returns -1, 0, or 1 depending on
whether p_ops precedes, is the same as, or follows q_ops.

\"Precede\" here means:
  - p precedes q (i.e., cmp_vec(p, q) returns -1), or
  - p equals q and length([x...]) is greater than length([y...])
  - p equals q, length([x...]) equals length([y...]), and [x...] precedes
    [y...] lexicographically.
"""
function cmp_pgroup(p_ops, q_ops)
    (p, u) = p_ops
    (q, v) = q_ops

    cmp_pq = cmp_vec(p, q)

    if (cmp_pq != 0)
        return cmp_pq
    end

    m = length(u)
    n = length(v)

    if (m != n)
        return (m > n) ? -1 : 1
    end

    return cmp_samelen(u, v)
end



"""
Define ordering of monomials, for printing and choosing canonical monomials
out of a group (e.g., using the minimum cycle to represent the trace).

If degree(x) < degree(y) then x precedes y. Otherwise, if x and y have the
same degree then the result is true if x precedes y lexicographically. The
way x is compared to y in that case favours larger groups of adjacent
operators sharing the same party. For example,

  A1 < A2 ,

but

  A2 A1 B1 < A1 B1 B2

because A2 A1 B1 starts with two As while A1 B1 B2 only starts with one. The
reason for this is to keep adjacent operators sharing the same party together
when using the lexicographically first cycle as the representative of the
trace. This way, for example, we represent trace(A2 A1 A_B1) with A2 A1 A_B1
rather than A1 A_B1 A2.
"""
function Base.isless(x::Monomial, y::Monomial)
    deg_x = degree(x)
    deg_y = degree(y)

    if (deg_x != deg_y)
        return (deg_x < deg_y)
    end

    if (deg_x == 0)
        return false
    end

    wx = x.word
    wy = y.word

    nx = length(wx)
    ny = length(wy)

    nx_le_ny = (nx < ny)

    m = (nx_le_ny ? nx : ny)
    j = 1

    while (j <= m)
        cmp_xy = cmp_pgroup(wx[j], wy[j])

        if (cmp_xy != 0)
            return (cmp_xy == -1)
        end

        j += 1
    end

    return nx_le_ny
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
given indices j and k.

parties_isect(p, q) returns the same thing as !isempty(intersect(p, q)), but
uses that the elements of p and q are strictly increasing to 
"""
function parties_isect(p::PartyVec, q::PartyVec)
    if ((m = length(p)) == 0) || ((n = length(q)) == 0)
        return false
    end

    (j, x) = (1, p[1])
    (k, y) = (1, q[1])

    while true
        if x == y
            return true
        elseif x < y
            if (j += 1) > m
                return false
            end
         
            x = p[j]
        else
            if (k += 1) > n
                return false
            end

            y = q[k]
        end
    end
end

"""
Test if two party vecs should be swapped (return :swap), merged
(return :join), or left in the order they already are (return :leave).
"""
function swap_or_join(p::PartyVec, q::PartyVec)
    cmp_pq = cmp_vec(p, q)

    if (cmp_pq == 1)
        return parties_isect(p, q) ? :leave : :swap
    end

    return (cmp_pq == 0) ? :join : :leave
end

"""
We have a vector y = [(q, [ops...])...] of groups of adjacent operators
[ops...] sharing the same parties q and 
"""
function insert_at(p::PartyVec, y::OpVector, n::Int)
    for k in 1:n
        q = y[k][1]
        result = swap_or_join(p, q)

        if result === :join
            return (k, :join)
        elseif result === :leave
            return (k, :leave)
        end
    end

    return (n+1, :leave)
end

"""
Re-minimise a word that has just had a new group of operators added to its
front. Input has to contain at least one element.

For example, suppose we had the monomial A1 B1 and we added A_C1 to the front
to get A_C1 A1 B1. remin! reorders the OpVector representing this to the
vector corresponding to B1 A_C1 A1, which is the rearrangement of A_c1 A1 B1
that comes first lexicographically.
"""
function remin!(word)
    n = length(word)

    p0 = first(word[1])
    parties = Set{Int}(p0)

    move_front = Int[]
    move_back = Int[1]

    for k in 2:n
        p = first(word[k])

        if parties_isect(parties, p) || (cmp_vec(p, p0) != -1)
            parties = union!(parties, p)
            push!(move_back, k)
        else
            push!(move_front, k)
        end
    end

    perm = vcat(move_front, move_back)
    permute!(word, perm)
    return
end

"""
Determine the product of two monomials, represented as lists
[(p, [ops...])...] of pairs of parties and operators.
"""
function join_words(x::OpVector, y::OpVector)
    if length(x) == 0
        return (1, y)
    elseif (n = length(y)) == 0
        return (1, x)
    end

    coeff = 1
    word = copy(y)
    left = copy(x)
    k = n

    while !isempty(left)
        pu = pop!(left)
        (p, u) = pu
        (k, action) = insert_at(p, word, n)

        if action === :join
            (c, w) = join_ops(u, word[k][2])

            if iszero(c)
                return (0, Id_word)
            end

            coeff *= c

            if !isempty(w)
                word[k] = (p, w)
                remin!(@view(word[k:end]))
            else
                append!(left, @view(word[1:(k-1)]))
                deleteat!(word, 1:k)
                n -= k
            end
        else
            insert!(word, k, pu)
            remin!(@view(word[k:end]))
            n += 1
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



function word_conj(word::OpVector)
    n = length(word)
    result = [(p, map(conj, Iterators.reverse(u)))
              for (p, u) in reverse(word)]

    for k in (n-1):-1:1
        remin!(@view(result[k:n]))
    end

    return result
end

Base.conj(m::Monomial) = Monomial(word_conj(m.word))
Base.adjoint(m::Monomial) = Monomial(word_conj(m.word))

function parties_isect(p::Set{Int}, q::PartyVec)
    for k in q
        if k in p
            return true
        end
    end

    return false
end

function isect_or_equal(p::PartyVec, q::PartyVec)
    if p == q
        return :equal
    elseif parties_isect(p, q)
        return :isect
    else
        return :disjoint
    end
end

function cjoin_or_leave(p::PartyVec, j::Int, word::OpVector, n::Int)
    k = n

    while k > j
        q = word[k][1]

        result = isect_or_equal(p, q)

        if result === :equal
            return (k, :join)
        elseif result === :isect
            return (k, :leave)
        end

        k -= 1
    end

    return (j, :join)
end

function min_party_cycle(word::OpVector, m=length(word))
    if m < 2
        return word
    end

    off_min = 0

    for off in 1:(m-1)
        if isless_samelen(word, word, off, off_min, m, m)
            off_min = off
        end
    end
    
    head = view(word, 1:off_min)
    tail = view(word, (off_min+1):m)

    return vcat(tail, head)
end

function ctrace(word::OpVector)
    coeff = 1
    word = copy(word)

    n = length(word)
    j = 1
    parties = Set{Int}()

    while j <= n
        (p, u) = word[j]

        if parties_isect(parties, p)
            j += 1
            continue
        end

        (k, result) = cjoin_or_leave(p, j, word, n)

        if result === :join
            if j == k
                (c, w) = trace(u)

                if c == 0
                    return (0, Id_word)
                end

                coeff *= c
                word[j] = (p, w)
                j += 1
                parties = union!(parties, p)
            else
                (c, w) = join_ops(word[k][2], u)

                if c == 0
                    return (0, Id_word)
                end

                coeff *= c
                deleteat!(word, k)
                n -= 1

                if isempty(w)
                    deleteat!(word, j)
                    n -= 1
                else
                    word[j] = (p, w)
                    j += 1
                    parties = union!(parties, p)
                end
            end
        else
            j += 1
            parties = union!(parties, p)
        end
    end

    return (coeff, min_party_cycle(word, n))
end

function ctrace(m::Monomial)
    (coeff, word) = ctrace(m.word)
    return (coeff, Monomial(word))
end
