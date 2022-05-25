Block = Dict{Integer}{SparseMatrixCSC}
Moments = Dict{Monomial}{Block}

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
Construct the NPA moment matrix.

The argument operators can in general be an array of arrays of operators
(blocks), e.g. [[Id], [A1 + B1], [A1 - B1]]. It can also be a simple array of
operators, in which case it is treated the same as an array containing a
single array of operators, e.g., [[Id, A1, A2]]). In either case the return
value is a dictionary with:

  * as keys: monomials obtained by multiplying operators in the same blocks
    together.

  * as values: dictionaries with integers (block numbers) as the keys and
    sparse matrices as the values. The sparse matrices contain coefficients
    obtained from the results of multiplying operators with the conjugates of
    operators in the same block together.

"""
function npa_moments(operators)
    moments = Moments()

    if isempty(operators)
        return moments
    end

    if first(operators) isa Union{Number,Monomial,Polynomial}
        operators = [operators]
    end

    for (nblock, ops) in enumerate(operators)
        N = length(ops)

        iops = collect(enumerate(ops))

        for (i, x) in iops
            for (j, y) in iops[i:end]
                p = Polynomial(conj_min(conj(x)*y))

                for (c, m) in p
                    if haskey(moments, m)
                        moment = moments[m]

                        if haskey(moment, nblock)
                            sparse_sym_set!(moment[nblock], i, j, c)
                        else
                            moment[nblock] = sparse_sym(N, i, j, c)
                        end
                    else
                        moments[m] = Block(nblock
                                           => sparse_sym(N, i, j, c))
                    end
                end
            end
        end
    end

    return moments
end

function Base.show(io::IO, moments::Moments)
    for (m, moment) in sort(moments)
        println(io)
        println(io, m, " => \n")

        for (b, mat) in sort(moment)
            println(io, "block ", b, ':')
            println(io, mat)
            println(io)
        end

        println(io)
    end
end

function Base.display(moments::Moments)
      for (m, moment) in sort(moments)
        println()
        println(m, " => \n")

        for (b, mat) in sort(moment)
            println("block ", b, ':')
            display(Array(mat))
            println()
        end
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
        moment0 = moments[m0]
        delete!(moments, m0)

        q = constraint[m0]
        constraint[m0] = 0

        for (c, m) in constraint
            moment = moments[m]

            for (b, g) in moment0
                delta = rdiv(c, q)*g

                if haskey(moment, b)
                    moment[b] -= delta
                else
                    moment[b] = -delta
                end
            end
        end
    end

    return (expr, moments)
end

function npa2sdp(expr,
                 constraints,
                 level)
    monomials = ops_at_level([expr, constraints], level)

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
    gamma = Dict()

    for (m, moment) in moments
        v = vars[m]

        for (b, g) in moment
            if haskey(gamma, b)
                gamma[b] += g*vars[m]
            else
                gamma[b] = g*vars[m]
            end
        end
    end

    constraints = [(g in :SDP) for g in values(gamma)]

    if goal in (:maximise, :maximize, :max)
        problem = maximize(objective, constraints)
    elseif goal in (:minimise, :minimize, :min)
        problem = minimize(objective, constraints)
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
