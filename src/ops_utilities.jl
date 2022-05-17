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

max_monomial(x) = maximum(all_monomials(x))



coefficients(p::Polynomial) = values(p.terms)



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

function substitute!(ps, x::Polynomial, m::Monomial)
    for p in ps
        substitute!(p, x, m)
    end

    return ps
end

function substitute!(ps::Set, x::Polynomial, m::Monomial)
    delete!(ps, x)

    for p in ps
        delete!(ps, p)
        p = substitute!(p, x, m)

        if !iszero(p)
            push!(ps, p)
        end
    end

    return ps
end

"Remove lexicographically highest monomial in x from p assuming <x> = 0."
function substitute!(p::Polynomial, x::Polynomial)
    return !iszero(x) ? substitute!(p, x, max_monomial(x)) : p
end

substitute(p::Polynomial, x, m) = substitute!(copy(p), x, m)
substitute(p::Polynomial, x) = substitute!(copy(p), x)


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


canonical(x::Number) = abs(x)
canonical(x::RNum) = 1
canonical(m::Monomial) = m

function canonical_factor(p::Polynomial)
    cfs = Rational[x for (_, x) in sort(p)
                       if (x isa RNum) && (!iszero(x))]

    if isempty(cfs)
        return 0
    end

    f = gcd(cfs)

    if !iszero(f) && (cfs[1] < 0)
        f = -f
    end

    return f
end

function canonical!(p::Polynomial)
    f = canonical_factor(p)

    if iszero(f)
        return p
    end

    return mul!(p, inv(f))
end

"Return polynomial scaled to make all its rational coefficients integers."
function canonical(p::Polynomial)
    f = canonical_factor(p)

    if iszero(f)
        return p
    end

    return inv(f)*p
end

function canonical_factor(p::Polynomial, m::Monomial)
    cfs = Rational[x for x in coefficients(p)
                       if (x isa RNum) && (!iszero(x))]

    if isempty(cfs)
        return 0
    end

    f = gcd(cfs)

    if !iszero(f) && (p[m] < 0)
        f = -f
    end

    return f
end

"""
Scale polynomial to make all its rational coefficients integers. Also make
p[m] positive if it is nonzero.
"""
function canonical!(p::Polynomial, m::Monomial)
    f = canonical_factor(p, m)

    if iszero(f)
        return p
    end

    return mul!(p, inv(f))
end

function max_constraint(constraints)
    (c0, rest) = Iterators.peel(constraints)
    m0 = max_monomial(c0)

    for c in constraints
        m = max_monomial(c)

        if m > m0
            m0 = m
            c0 = c
        end
    end

    return (m0, c0)
end

Linspace = Dict{Monomial,Polynomial}

"""
Add polynomial to space if it is not a linear combination of the polynomials
already in the space, assuming space is already in a certain reduced form.
space is a table of polynomials together with the lexicographically maximal
monomial in each polynomial.
Reduced form means that the highest monomial m in each polynomial space[m] in
space appears in no other polynomial (i.e., has been eliminated) in space.
Specifically,
  space[m][m] is nonzero,
and
  space[m][m1] = 0
for m, m1 in keys(space).
"""
function extend!!(space::Linspace, p::Polynomial)
    if iszero(p)
        return space
    end

    # Remove highest monomials already in space from p.
    # If p becomes 0 this means p is not linearly independent of space and
    # we can stop.

    for (m, q) in space
        if iszero(substitute!(p, q, m))
            return space
        end
    end

    # Now remove highest remaining monomial m0 of p from polynomials in
    # space.

    m0 = max_monomial(p)
    p = canonical!(p, m0)

    for (m, q) in space
        if !iszero(q[m0])
            canonical!(substitute!(q, p, m0), m)
        end
    end

    space[m0] = p

    return space
end

function extend!(space::Linspace, p::Polynomial)
    return extend!!(space, copy(p))
end



function linspace(operators)
    space = Linspace()

    for o in operators
        extend!(space, Polynomial(o))
    end

    return space
end



function reduce_expr!(expr, space::Linspace)
    return reduce_expr!(Polynomial(expr), space)
end

function reduce_expr!(expr::Polynomial, space::Linspace)
    for (m, p) in space
        substitute!(expr, p, m)
    end

    return expr
end



function reduce_expr(expr, space::Linspace)
    return reduce_expr!(new_polynomial(expr), space)
end
