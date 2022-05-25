Moments = Dict{Monomial}{SparseMatrixCSC}

function sparse_sym_set!(matrix, i, j, val)
    matrix[i, j] = val

    if i != j
        matrix[j, i] = val
    end
end

function sparse_sym(N, i, j, val)
    if i == j
        return sparse([i], [i], [val], N, N)
    else
        return sparse([i, j], [j, i], [val, val], N, N)
    end
end

"""
Construct the NPA moment matrix. Returns a dictionary with monomials
appearing in the moment matrix as keys and sparse matrices with 1s at
their positions.
"""
function npa_moments(operators)
    moments = Moments()

    N = length(operators)

    ops = collect(enumerate(operators))

    for (i, x) in ops
        for (j, y) in ops[i:end]
            p = Polynomial(conj_min(conj(x)*y))

            for (c, m) in p
                if haskey(moments, m)
                    sparse_sym_set!(moments[m], i, j, c)
                else
                    moments[m] = sparse_sym(N, i, j, c)
                end
            end
        end
    end

    return moments
end

function Base.display(moments::Moments)
    for (m, moment) in moments
        println()
        println(m, " => ")
        display(Array(moment))
    end
end



default_solver = SCS.Optimizer

function set_solver(solver)
    global default_solver = solver
end



"""
Generate the NPA relaxation for a given quantum optimisation problem (an
operator expr whose expectation we want to maximise with the expectation
values of the operators constraints set to zero).
"""
function npa2sdp(expr,
                 constraints,
                 moments::Moments)
    # Reduce constraints to canonical form

    if !(constraints isa Linspace)
        constraints = linspace(constraints)
    end

    constraints = Dict(m => conj_min(p) for (m, p) in constraints)

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

        for (c, m) in constraint
            moments[m] -= rdiv(c, q)*G
        end
    end

    return (expr, moments)
end

function npa2sdp(expr,
                 constraints,
                 level)
    monomials = ops_at_level(level, [expr, constraints])

    return npa2sdp(expr,
                   constraints,
                   npa_moments(monomials))
end

"""
Convert an SDP returned by npa2sdp to the Convex.jl problem format.
"""
function sdp2Convex(expr, moments; goal=:maximise)
    vars = Dict(m => ((m == Id) ? 1 : Variable())
                for m in keys(moments))

    objective = sum(c*vars[m] for (c, m) in expr)
    gamma = sum(g*vars[m] for (m, g) in moments)

    if goal in (:maximise, :maximize, :max)
        problem = maximize(objective, [(gamma in :SDP)])
    elseif goal in (:minimise, :minimize, :min)
        problem = minimize(objective, [(gamma in :SDP)])
    end

    return problem
end

function npa_opt(expr,
                 constraints,
                 level_or_moments;
                 solver=default_solver,
                 verbose=true,
                 goal=:maximise)
    (expr, moments) = npa2sdp(expr, constraints, level_or_moments)

    problem = sdp2Convex(expr, moments, goal=goal)
    solve!(problem, solver, verbose=verbose, silent_solver=true)

    return problem.optval
end



function npa_max(expr, level; solver=default_solver, verbose=true)
    return npa_opt(expr, [], level,
                   solver=solver,
                   verbose=verbose,
                   goal=:maximise)
end

function npa_max(expr, constraints, level;
                 solver=default_solver,
                 verbose=true)
    return npa_opt(expr, constraints, level,
                   solver=solver,
                   verbose=verbose,
                   goal=:maximise)
end

function npa_min(expr, level; solver=default_solver, verbose=true)
    return npa_opt(expr, [], level,
                   solver=solver,
                   verbose=verbose,
                   goal=:minimise)
end

function npa_min(expr, constraints, level;
                 solver=default_solver,
                 verbose=true)
    return npa_opt(expr, constraints, level,
                   solver=solver,
                   verbose=verbose,
                   goal=:minimise)
end
