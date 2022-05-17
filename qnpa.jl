using Base.Iterators
#  flatten, zip (filter?)

using Combinatorics
#  powerset

using Convex
using SCS

using SparseArrays



include("src/operators.jl")
include("src/ops_predefined.jl")



# NPA

Moments = Dict{Monomial}{SparseMatrixCSC}

function npa_moments(monomials)
    moments = Moments()

    N = length(monomials)

    ops = collect(enumerate(monomials))

    for (i, x) in ops
        for (j, y) in ops[i:end]
            m = conj(x)*y

            if m == 0
                continue
            end

            m = min(m, conj(m))

            if haskey(moments, m)
                moments[m][i, j] = 1

                if i != j
                    moments[m][j, i] = 1
                end
            else
                if i == j
                    moments[m] = sparse([i], [i], [1], N, N)
                else
                    moments[m] = sparse([i, j], [j, i], [1, 1], N, N)
                end
            end
        end
    end

    return moments
end


function ops_at_level(n::Integer, source)
    ops = operators(source)
    return ops_at_level(n, ops)
end

function ops_at_level(n::Integer, ops::Set{Monomial})
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

function ops_at_level(n::AbstractString, source)
    ops = operators(source, by_party=true)
    return ops_at_level(n, ops)
end

function ops_at_level(n::AbstractString, ops::Dict{Integer,Set{Monomial}})
    result = Set{Monomial}()

    for term in split(n, '+')
        ms = make_term(term, ops)
        result = union!(result, ms)
    end

    return sort(result)
end



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

function nonzero_products(s1, s2)
    return filter(!iszero, Set(x*y for x in s1 for y in s2))
end

function monomial_products(monomials, power::Int)
    @assert power > 0

    result = Set([Id])

    while power > 0
        result = nonzero_products(result, monomials)
        power -= 1
    end

    return result
end

function make_term(term, ops::Dict{Integer,Set{Monomial}})
    level = parse_level_term(term)

    if level isa Integer
        all_ops = union(values(ops)...)

        return ops_at_level(level, all_ops)
    else
        result = Set([Id])
        
        for (party, power) in level
            ms = monomial_products(ops[party], power)
            result = nonzero_products(result, ms)
        end

        return result
    end
end



default_solver = SCS.Optimizer

function set_solver(solver)
    global default_solver = solver
end

function npa_opt(expr,
                 constraints,
                 moments::Moments;
                 solver=default_solver,
                 goal=:maximise)
    # Reduce constraints to canonical form

    if !(constraints isa Linspace)
        constraints = linspace(constraints)
    else
        constraints = deepcopy(constraints)
    end

    if haskey(constraints, Id)
        @error "Contradiction Id = 0 in constraints."
    end

    # Reduce the objective expression, using constraints to eliminate
    # monomials
    expr = reduce_expr(expr, constraints)

    moments = deepcopy(moments)

    for (m0, constraint) in constraints
        G = moments[m0]
        delete!(moments, m0)

        q = constraint[m0]
        constraint[m0] = 0

        for (m, c) in constraint
            moments[m] -= rdiv(c, q)*G
        end
    end

    vars = Dict(m => ((m == Id) ? 1 : Variable())
                for m in keys(moments))

    objective = sum(c*vars[m] for (m, c) in expr)
    gamma = sum(g*vars[m] for (m, g) in moments)

    if goal in (:maximise, :maximize, :max)
        problem = maximize(objective, [(gamma in :SDP)])
    elseif goal in (:minimise, :minimize, :min)
        problem = minimize(objective, [(gamma in :SDP)])
    end

    solve!(problem, solver, silent_solver=true)

    return evaluate(objective)
end

function npa_opt(expr,
                 constraints,
                 level;
                 solver=default_solver,
                 goal=:maximise)
    monomials = ops_at_level(level, [expr, constraints])

    return npa_opt(expr,
                   constraints,
                   npa_moments(monomials),
                   solver=solver,
                   goal=goal)
end



function npa_max(expr, level; solver=default_solver)
    return npa_opt(expr, [], level, solver=solver, goal=:maximise)
end

function npa_max(expr, constraints, level; solver=default_solver)
    return npa_opt(expr, constraints, level,
                   solver=solver,
                   goal=:maximise)
end

function npa_min(expr, level; solver=default_solver)
    return npa_opt(expr, [], level, solver=solver, goal=:minimise)
end

function npa_min(expr, constraints, level; solver=default_solver)
    return npa_opt(expr, constraints, level,
                   solver=solver,
                   goal=:minimise)
end
