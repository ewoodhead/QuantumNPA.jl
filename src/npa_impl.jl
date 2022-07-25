Moments = Dict{Monomial}{BlockDiagonal}

function sparse_sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end
end

function sparse_sym(N, i, j, val)
    if i == j
        return sparse([i], [i], [val], N, N)
    else
        return sparse([i, j], [j, i], [val, val], N, N)
    end
end

function npa_moments_block(operators; f=identity)
    N = length(operators)
    iops = collect(enumerate(operators))
    block = Dict{Monomial,SparseMatrixCSC}()

    for (i, x) in iops
        for (j, y) in iops[i:end]
            z = conj_min(conj(x)*y, f=f)
            p = Polynomial(z)

            for (c, m) in p
                if !haskey(block, m)
                    block[m] = sparse_sym(N, i, j, c)
                else
                    sparse_sym_add!(block[m], i, j, c)
                end
            end
        end
    end

    return block
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
function npa_moments(operators; f=identity)
    if isempty(operators)
        return moments
    end

    if first(operators) isa Union{Number,Monomial,Polynomial}
        operators = [operators]
    end

    nblocks = length(operators)
    bsizes = length.(operators)
    blocks = [npa_moments_block(o, f=f) for o in operators]

    ms = monomials(keys(block) for block in blocks)

    moments = Moments()

    for m in ms
        blocks_m = [(haskey(block, m)
                     ? block[m]
                     : (n -> spzeros(n, n))(bsizes[b]))
                    for (b, block) in enumerate(blocks)]

        moments[m] = BlockDiagonal(blocks_m)
    end

    return moments
end



function SparseArrays.dropzeros!(matrix::BlockDiagonal)
    for blk in blocks(matrix)
        dropzeros!(blk)
    end

    return matrix
end

"""
Generate the NPA relaxation for a given quantum optimisation problem (an
operator expr whose expectation we want to maximise with the expectation
values of the operators constraints set to zero).
"""
function npa2sdp(expr,
                 constraints,
                 moments::Moments;
                 f=identity)
    # Reduce constraints to canonical form

    expr = conj_min(expr; f=f)
    constraints = linspace(map(m -> conj_min(m, f=f), constraints))

    if haskey(constraints, Id)
        @error "Contradiction Id = 0 in constraints."
    end

    # Reduce the objective expression, using constraints to eliminate
    # monomials
    expr = reduce_expr(expr, constraints)
    moments = deepcopy(moments)

    for (m0, constraint) in constraints
        q = constraint[m0]
        constraint[m0] = 0

        moment0 = moments[m0]
        delete!(moments, m0)

        for (c, m) in constraint
            moments[m] -= rdiv(c, q) * moment0
        end
    end

    # Remove any zero coefficients that might be stored explicitly in the
    # sparse matrices in the blocks.
    # for matrix in values(moments)
    #    dropzeros!(matrix)
    # end
    
    moments = Moments(m => mat
                      for (m, mat) in moments
                          if !iszero(mat))

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

npa2sdp(expr, lvl_or_constraints) = npa2sdp(expr, [], lvl_or_constraints)



function SparseArrays.findnz(mat::BlockDiagonal)
    base = 0

    is = Int[]
    js = Int[]
    cs = []

    for blk in blocks(mat)
        (is1, js1, cs1) = findnz(blk)
        append!(is, is1 .+ base)
        append!(js, js1 .+ base)
        append!(cs, cs1)
        base += first(size(mat))
    end

    return (is, js, cs)
end

function firstnz(moment::BlockDiagonal)
    base = 0

    for blk in blocks(moment)
        if !iszero(blk)
            (i, j, c) = first((i, j, c)
                              for (i, j, c) in zip(findnz(blk)...))
            @assert !iszero(c)

            if (i > j)
                (i, j) = (j, i)
            end
            
            return (i + base, j + base, c)
        else
            base += first(size(blk))
        end
    end
end


function bspzeros(bsizes)
    return BlockDiagonal([spzeros(n, n) for n in bsizes])
end

function Base.zero(bm::BlockDiagonal)
    return bspzeros(first.(blocksizes(bm)))
end

function BlockDiagonals.blocksizes(moments::Moments)
    if isempty(moments)
        return []
    else
        return first.(blocksizes(first(moments)[2]))
    end
end



if !@isdefined(default_solver)
    default_solver = SCS.Optimizer
end

function set_solver(solver)
    global default_solver = solver
end

function set_verbosity!(model, verbose)
    if !isnothing(verbose)
        (!verbose ? set_silent : unset_silent)(model)
    end
end



function sdp2jump(expr, moments;
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing,
                  fixed=Id)
    if goal in (:maximise, :maximize, :max)
        maximise = true
        s = 1
    elseif goal in (:minimise, :minimize, :min)
        maximise = false
        s = -1
    end
    
    model = !isnothing(solver) ? Model(solver) : Model()

    Z = [@variable(model, [1:n, 1:n], PSD) for n in blocksizes(moments)]

    objective = (sum(LinearAlgebra.tr(s*G*Z[b])
                     for (b, G) in enumerate(blocks(moments[fixed])))
                 + expr[fixed])
    
    if maximise
        @objective(model, Min, objective)
    else
        @objective(model, Max, objective)
    end

    for (m, moment) in moments
        if m != fixed
            c = expr[m]
            
            @constraint(model,
                        sum(LinearAlgebra.tr(F*Z[b])
                            for (b, F) in enumerate(blocks(moment)))
                        + s*c == 0)
        end
    end

    set_verbosity!(model, verbose)

    return model
end



function npa2jump(expr,
                  constraints,
                  level_or_moments;
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing)
    (expr, moments) = npa2sdp(expr, constraints, level_or_moments)
    model = sdp2jump(expr, moments,
                     goal=goal,
                     solver=solver,
                     verbose=verbose)

    return model
end

function npa2jump(expr, level_or_moments;
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing)
    return npa2jump(expr, [], level_or_moments,
                    goal=goal,
                    solver=solver,
                    verbose=verbose)
end



function npa_opt(expr,
                 constraints,
                 level_or_moments;
                 goal=:maximise,
                 solver=default_solver,
                 verbose=false)
    model = npa2jump(expr, constraints, level_or_moments, goal=goal)

    set_optimizer(model, solver)

    if !verbose
        set_silent(model)
    end

    optimize!(model)

    return objective_value(model)
end



function npa_max(expr, constraints, level;
                 solver=default_solver,
                 verbose=false)
    return npa_opt(expr, constraints, level,
                   goal=:maximise,
                   solver=solver,
                   verbose=verbose)
end

function npa_max(expr, level; solver=default_solver, verbose=false)
    return npa_max(expr, [], level,
                   solver=solver,
                   verbose=verbose)
end



function npa_min(expr, constraints, level;
                 solver=default_solver,
                 verbose=false)
    return npa_opt(expr, constraints, level,
                   goal=:minimise,
                   solver=solver,
                   verbose=verbose)
end

function npa_min(expr, level; solver=default_solver, verbose=false)
    return npa_min(expr, [], level,
                   solver=solver,
                   verbose=verbose)
end
