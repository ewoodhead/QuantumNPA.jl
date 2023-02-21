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



# Implement lexicographical comparison of party vectors and monomials.

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

isless_pgroup(p_ops, q_ops) = (cmp_pgroup(p_ops, q_ops) == -1)

function isless_word(wx::OpVector, wy::OpVector)
    nx = length(wx)
    ny = length(wy)

    nx_le_ny = (nx < ny)

    m = (nx_le_ny ? nx : ny)

    for j in 1:m
        cmp_xy = cmp_pgroup(wx[j], wy[j])

        if (cmp_xy != 0)
            return (cmp_xy == -1)
        end
    end

    return nx_le_ny
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

    return isless_word(x.word, y.word)
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
    return word
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



# Implementation of conjugate/adjoint of a monomial.

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



# Helper functions to reduce monomials under trace. The actual implementation
# of trace() is deferred to the file ops_polynomial.jl since this code
# supports the possibility that the trace might have to be a polynomial,
# specifically, a monomial multiplied by a coefficient different from 1.

"""
Test if p and q intersect. Assumes that q is a vector of integers in stricty
increasing order.
"""
function parties_isect(p::Set{Int}, q::PartyVec)
    for k in q
        if k in p
            return true
        end
    end

    return false
end



"""
Find if a group (q, [...]) of operators in a word = [(p, [ops...])...] can
be brought to the front then merged with an operator in the same word going
back through it from the end. For example, in the word corresponding to
A1 A_B1 A2, this would determine that the A1 at the beginning could be merged
with the A2 at the end.

Parameters:

p : party vector at position j
j : position of operator group in word
word : the entire word
n : the length of word

Returns a pair (k, bool) indicating whether and where the operator group
can be joined. Testing for j == k means the group can be joined with itself
(so something like A1 A2 A1 can be simplified to A2).
"""
function can_join_back(p::PartyVec, j::Int, word::OpVector, n::Int)
    for k in n : -1 : (j+1)
        q = word[k][1]

        if p == q
            return (k, true)
        elseif parties_isect(p, q)
            return (k, false)
        end
    end

    return (j, true)
end

"""
Join operator group at index j in word with a group taken from later in word,
if possible.
"""
function join_back!(j::Int, parties::Set{Int}, word::OpVector, n::Int)
    (p, u) = word[j]

    if parties_isect(parties, p)
        return (j+1, n, 1, word)
    end

    (k, result) = can_join_back(p, j, word, n)

    if !result
        union!(parties, p)
        return (j+1, n, 1, word)
    end

    (coeff, w) = ((j == k) ? trace(u) : join_ops(word[k][2], u))

    if coeff == 0
        return (1, 0, 0, Id_word)
    end

    if (j != k)
        deleteat!(word, k)
        n -= 1

        if isempty(w)
            deleteat!(word, j)
            return (j, n-1, coeff, word)
        end
    end

    word[j] = (p, w)
    union!(parties, p)

    return (j+1, n, coeff, word)
end

"""
Join operator groups from the back of word to operator groups from the front,
where possible.

For example, given the word corresponding to A1 B1 A_B1 B2, this function
returns the word corresponding to A1 B2 B1 A_B1.
"""
function join_back(word::OpVector)
    coeff = 1
    word = copy(word)
    n = length(word)
    parties = Set{Int}()
    j = 1

    while j <= n
        (j, n, c, word) = join_back!(j, parties, word, n)
        coeff *= c
    end

    return (coeff, word)
end

"Return minimal cycle of groups of operators in word."
function min_cycle(word::OpVector)
    n = length(word)
    n1 = n - 1
    mcycle = word

    for _ in 1:n1
        word = remin!(vcat(word[end], word[1:n1]))

        if isless_word(word, mcycle)
            mcycle = word
        end
    end

    return mcycle
end

function ctrace(word::OpVector)
    (coeff, word) = join_back(word)
    return (coeff, min_cycle(word))
end

function ctrace(m::Monomial)
    (coeff, word) = ctrace(m.word)
    return (coeff, Monomial(word))
end
