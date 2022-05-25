"""
Return (without duplications) all nonzero products of elements in s1 and s2.
"""
function nonzero_products(s1, s2)
    return filter(!iszero, Set(x*y for x in s1 for y in s2))
end

"Return all unique nonzero n-fold products of given monomials."
function monomial_products(monomials, power::Int)
    @assert power > 0

    result = Set([Id])

    while power > 0
        result = nonzero_products(result, monomials)
        power -= 1
    end

    return result
end

"""
Parse and return representation of a term contributing to an NPA hierarchy
level. This returns an integer if the term contains only digits and a
dictionary of party numbers and the powers to which they appear otherwise.

For example:

  parse_level_term("2")           =>  2

  parse_level_term("A^2 B B E")   =>  Dict(1 => 2,
                                           2 => 2,
                                           5 => 1)
"""
function parse_level_term(term)
    term = strip(term)

    if all(isdigit, term)
        return parse(Int, term)
    end

    contribs = Dict{Int,Int}()

    for party_str in split(term)
        if '^' in party_str
            (party_str, p_str) = split(party_str, '^')
            power = parse(Int, p_str)
        else
            power = 1
        end

        party = party_num(party_str)

        if haskey(contribs, party)
            contribs[party] += power
        else
            contribs[party] = power
        end
    end

    return contribs
end

"""
Return a set of operators corresponding to a particular level term, such as
"1" or "A B".
"""
function make_term(term, ops::Dict{Integer,Set{Monomial}})
    level = parse_level_term(term)

    if level isa Integer
        all_ops = union(values(ops)...)

        return ops_at_level(all_ops, level)
    else
        result = Set([Id])
        
        for (party, power) in level
            ms = monomial_products(ops[party], power)
            result = nonzero_products(result, ms)
        end

        return result
    end
end



"Return operators taken fron source at a given level of the NPA hierarchy."
function ops_at_level(source, n::Integer)
    ops = operators(source)
    return ops_at_level(ops, n)
end

function ops_at_level(ops::Set{Monomial}, n::Integer)
    @assert n >= 0

    ops = copy(ops)
    push!(ops, Id)
    result = Set([Id])

    while n > 0
        result = filter(!iszero, Set(x*y for x in result for y in ops))
        n -= 1
    end

    return sort(result)
end

function ops_at_level(source, n::AbstractString)
    ops = operators(source, by_party=true)
    return ops_at_level(ops, n)
end

function ops_at_level(ops::Dict{Integer,Set{Monomial}}, n::AbstractString)
    result = Set{Monomial}()

    for term in split(n, '+')
        ms = make_term(term, ops)
        result = union!(result, ms)
    end

    return sort(result)
end
