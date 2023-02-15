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



function isless_samelen(x::OpVector,
                        y::OpVector,
                        offset_x::Int=0,
                        offset_y::Int=0,
                        lenx::Int=length(x),
                        leny::Int=length(y))
    if iszero(lenx)
        return false
    end

    j = 0
    jc = 1 + (j + offset_x) % lenx
    (p, u) = x[jc]

    k = 0
    kc = 1 + (k + offset_y) % leny
    (q, v) = y[kc]

    if (p != q)
        return (p < q)
    end

    mu = length(u)
    iu = 1
    ou = u[iu]

    mv = length(v)
    iv = 1
    ov = v[iv]

    while true
        if ou != ov
            return ou < ov
        end

        if iu < mu
            if iv < mv
                iu += 1
                ou = u[iu]
                iv += 1
                ov = v[iv]
            else
                k += 1
                kc = 1 + (k + offset_y) % leny
                q = y[kc][1]
                return p < q
            end
        else
            j += 1
            jc = 1 + (j + offset_x) % lenx

            if iv < mv
                p = x[jc][1]
                return p < q
            elseif j < len
                (p, u) = x[jc]

                k += 1
                kc = 1 + (k + offset_y) % leny
                (q, v) = y[kc]

                if p != q
                    return p < q
                end

                mu = length(u)
                iu = 1
                ou = u[iu]

                mv = length(v)
                iv = 1
                ov = v[iv]
            else
                return false
            end
        end
    end
end

function Base.isless(x::Monomial, y::Monomial)
    ox, oy = degree(x), degree(y)

    if ox != oy
        return ox < oy
    end

    return isless_samelen(x.word, y.word)
end



Base.zero(m::Monomial) = 0



conj_min(x::Number; f=identity) = real(f(x))

function conj_min(m::Monomial; f=identity)
    z = f(m)

    if iszero(z)
        return 0
    end

    return min(z, f(conj(m)))



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

function join_words(x::OpVector, y::OpVector)
    if (m = length(x)) == 0
        return (1, y)
    elseif (n = length(y)) == 0
        return (1, x)
    end

    coeff = 1
    word = copy(y)
    k = n

    for (j, pu) in Iterators.reverse(enumerate(x))
        (p, u) = pu
        (k, action) = insert_at(p, word, k)

        if action === :join
            (c, w) = join_ops(u, word[k][2])

            if iszero(c)
                return (0, Id_word)
            end

            coeff *= c

            if !isempty(w)
                word[k] = (p, w)
                k -= 1
            else
                deleteat!(word, k)
                n -= 1
                k = n
            end
        else
            insert!(word, k, pu)
            k -= 1
            n += 1
        end
    end

    return (coeff, word)
end

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

function join_words(x::OpVector, y::OpVector)
    if (m = length(x)) == 0
        return (1, y)
    elseif (n = length(y)) == 0
        return (1, x)
    end

    coeff = 1
    word = copy(y)
    k = n

    for (j, pu) in Iterators.reverse(enumerate(x))
        (p, u) = pu
        (k, action) = insert_at(p, word, k)

        if action === :join
            (c, w) = join_ops(u, word[k][2])

            if iszero(c)
                return (0, Id_word)
            end

            coeff *= c

            if !isempty(w)
                word[k] = (p, w)
                k -= 1
            else
                deleteat!(word, k)
                n -= 1
                k = n
            end
        else
            insert!(word, k, pu)
            k -= 1
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
    result = OpVector()

    for (p, u) in word
        w = map(conj, Iterators.reverse(u))

        k = 1
        n = length(result)

        while k <= n
            q = result[k][1]

            if swap_or_join(p, q) !== :swap
                break
            end
            
            k += 1
        end

        insert!(result, k, (p, w))
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
