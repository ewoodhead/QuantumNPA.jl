function sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end

    return matrix
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

  * as values: block-diagonal sparse matrices with coefficients obtained
    from multiplying the input operators together.
"""
function npa_moment(operators::Vector{<:Union{Monomial,Polynomial}})
    N = length(operators)
    iops = collect(enumerate(operators))
    moment = Polynomial((N, N))

    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(real_rep(conj(x)*y))

            for (c, m) in p
                if !hasmonomial(moment, m)
                    moment[m] = sym_add!(spzeros(N, N), i, j, c)
                else
                    sym_add!(moment[m], i, j, c)
                end
            end
        end
    end

    return moment
end

#function npa_moment(operators::Vector{Vector{T}} where T)
#    moment = npa_moment.(operators)
#    return blockdiag(moment, (sz) -> spzeros(Float64, sz))
#end

npa_moment(source, level) = npa_moment(npa_level(source, level))



"""
Generate the NPA relaxation for a given quantum optimisation problem (an
operator expr whose expectation we want to maximise with the expectation
values of the operator constraints set to zero).

The result is
"""
function npa2sdp(expr, level; eq=[], ge=[])
    moment = npa_moment([expr, eq, ge], level)
    return npa2sdp(expr, moment, eq=eq, ge=ge)
end

function npa2sdp(expr, moment::Polynomial; eq=[], ge=[])
    # Reduce constraints to canonical form
    expr = real_rep(expr)
    eq = linspace(map(real_rep, eq))
    ge = map(real_rep, ge)

    if haskey(eq, Id)
        @error "Contradiction Id = 0 in equality constraints."
    end

    # Reduce the objective expression, using constraints to eliminate
    # monomials.
    expr = reduce_expr(expr, eq)

    # Reduce moments using equality constraints.
    moment = reduce_expr(moment, eq)

    # Reduce inequality constraints then include them as inequalities along
    # with the original moment matrix.
    ge = reduce_exprs(ge, eq)

    return (expr, vcat([moment], ge))
end



if !@isdefined(default_solver)
    default_solver = SCS.Optimizer
end

"Set the default solver called by QuantumNPA"
function set_solver!(solver)
    global default_solver = solver
end



if !@isdefined(default_silent)
    default_silent = true
end

"Set default JuMP model verbosity"
function set_verbosity!(silent::Bool)
    global default_silent = silent
end

"Set JuMP model verbosity"
function set_verbosity!(model::Model, silent=default_verbosity)
    if !isnothing(silent)
        (silent ? set_silent : unset_silent)(model)
    end
end



function jump_model(solver=default_solver, silent=default_silent)
    model = (!isnothing(solver) ? Model(solver) : Model())

    if !isnothing(silent)
        set_verbosity!(model, silent)
    end

    return model
end



# This is required so that dot used in sdp2jump_d has acceptable performance
import LinearAlgebra.dot

function dot(A::SparseMatrixCSC, B::Symmetric{<:JuMP._MA.AbstractMutable})
    acc = zero(eltype(B))
    for j in 1:size(A, 2)
        for k in nzrange(A, j)
            add_to_expression!(acc, nonzeros(A)[k], B[rowvals(A)[k], j])
        end
    end
    return acc
end

dot(A::Symmetric{<:JuMP._MA.AbstractMutable}, B::SparseMatrixCSC) = dot(B,A)



function sdp2jump_d(expr, ineqs;
                    goal=:maximise,
                    solver=default_solver,
                    silent=default_silent)
    if goal in (:maximise, :maximize, :max)
        maximise = true
        s = 1
    elseif goal in (:minimise, :minimize, :min)
        maximise = false
        s = -1
    end

    model = jump_model(solver, silent)

    Zs = [@variable(model, [1:m, 1:n], PSD)
          for (m, n) in size_as_pair.(ineqs)]

    Ids = (ineq[Id] for ineq in ineqs)
    objective = (s*sum(dot(m,z)
                     for (m, z) in zip(Ids, Zs))
                 + expr[Id])

    if maximise
        @objective(model, Min, objective)
    else
        @objective(model, Max, objective)
    end

    mons = collect(m for m in monomials(expr, ineqs) if !isidentity(m))

    for m in mons
        c = expr[m]
        Fs = (ineq[m] for ineq in ineqs)
        tr_term = sum(dot(F,Z)
                      for (F, Z) in zip(Fs, Zs))

        @constraint(model, tr_term + s*c == 0)
    end

    return model
end



function npa2jump(expr, level_or_moments; eq=[], ge=[], kw...)
    (expr, moments) = npa2sdp(expr, level_or_moments, eq=eq, ge=ge)
    model = sdp2jump_d(expr, moments; kw...)
    return model
end



function npa_opt(expr, level_or_moments; eq=[], ge=[], kw...)
    model = npa2jump(expr, level_or_moments, eq=eq, ge=ge; kw...)
    optimize!(model)
    return objective_value(model)
end



npa_max(expr, level; kw...) = npa_opt(expr, level; goal=:maximise, kw...)
npa_min(expr, level; kw...) = npa_opt(expr, level; goal=:minimise, kw...)
