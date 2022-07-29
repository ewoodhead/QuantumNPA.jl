Moments = Dict{Monomial}{BlockDiagonal}

function sparse_sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end
end

function sparse_Nt(I, J, V, M, N)
    return SparseMatrixCSC{Float64}{Int}(sparse(I, J, V, M, N))
end

function sparse_sym(N, i, j, val)
    if i == j
        return sparse_Nt([i], [i], [val], N, N)
    else
        return sparse_Nt([i, j], [j, i], [val, val], N, N)
    end
end

function npa_moments_block(operators)
    N = length(operators)
    iops = collect(enumerate(operators))
    block = Dict{Monomial,SparseMatrixCSC}()

    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(conj_min(conj(x)*y))

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
function npa_moments(operators)
    if isempty(operators)
        return moments
    end

    if first(operators) isa Union{Number,Monomial,Polynomial}
        operators = [operators]
    end

    nblocks = length(operators)
    bsizes = length.(operators)
    blocks = npa_moments_block.(operators)

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

sp1x1(x) = (iszero(x)
            ? spzeros(Float64, 1, 1)
            : SparseMatrixCSC{Float64,Int}(sparse([1], [1], x)))

"""
Generate the NPA relaxation for a given quantum optimisation problem (an
operator expr whose expectation we want to maximise with the expectation
values of the operators constraints set to zero).
"""
function npa2sdp(expr,
                 level_or_moments;
                 eq=[],
                 ge=[])
    if level_or_moments isa Moments
        moments = level_or_moments
    else
        moments = npa_moments(ops_at_level([expr, eq, ge],
                                           level_or_moments))
    end
    
    # Reduce constraints to canonical form
    expr = conj_min(expr)
    eq = linspace(map(conj_min, eq))
    ge = map(conj_min, ge)

    if haskey(eq, Id)
        @error "Contradiction Id = 0 in equality constraints."
    end

    # Reduce the objective expression, using constraints to eliminate
    # monomials
    expr = reduce_expr(expr, eq)
    moments = deepcopy(moments)

    # Reduce moments using equality constraints.
    for (m0, constraint) in eq
        constraint = copy(constraint)
        q = constraint[m0]
        constraint[m0] = 0

        moment0 = moments[m0]
        delete!(moments, m0)

        for (c, m) in constraint
            moments[m] -= rdiv(c, q) * moment0
        end
    end

    # Reduce inequality constraints then absorb them into the moment matrix.
    # Basically, take the coefficients in the inequalities and add them as
    # 1x1 blocks to the moments.
    ge = [reduce_expr(ineq, eq) for ineq in ge]

    for (m, moment) in moments
        append!(blocks(moment),
                [sp1x1(ineq[m]) for ineq in ge])
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

function set_solver!(solver)
    global default_solver = solver
end

function set_verbosity!(model, verbose)
    if !isnothing(verbose)
        (!verbose ? set_silent : unset_silent)(model)
    end
end



function expr2objective(expr, vars)
    return expr[Id] + sum(c*vars[m] for (c, m) in expr if m != Id)
end

"""
Convert moments returned by npa2sdp() to moments in a format used by JuMP.jl
or Convex.jl.
"""
function moments2gamma(moments, vars)
    if isempty(moments)
        return []
    end

    n = nblocks(first(moments)[2])
    gamma = Vector(undef, n)

    for (m, moment) in moments
        var = ((m != Id) ? vars[m] : 1)

        for (b, g) in enumerate(blocks(moment))
            if isassigned(gamma, b)
                gamma[b] += g*var
            else
                gamma[b] = g*var
            end
        end
    end

    return gamma
end



function sdp2jump(expr, moments;
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing)
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
                     for (b, G) in enumerate(blocks(moments[Id])))
                 + expr[Id])
    
    if maximise
        @objective(model, Min, objective)
    else
        @objective(model, Max, objective)
    end

    for (m, moment) in moments
        if m != Id
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



function npa2jump(expr, level_or_moments;
                  eq=[],
                  ge=[],
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing)
    (expr, moments) = npa2sdp(expr, level_or_moments, eq=eq, ge=ge)

    model = sdp2jump(expr, moments,
                     goal=goal,
                     solver=solver,
                     verbose=verbose)

    return model
end



function npa_opt(expr, level_or_moments;
                 eq=[],
                 ge=[],
                 goal=:maximise,
                 solver=default_solver,
                 verbose=false)
    model = npa2jump(expr, level_or_moments,
                     eq=eq,
                     ge=ge,
                     goal=goal)

    set_optimizer(model, solver)

    if !verbose
        set_silent(model)
    end

    optimize!(model)

    return objective_value(model)
end



npa_max(expr, level; kw...) = npa_opt(expr, level; goal=:maximise, kw...)
npa_min(expr, level; kw...) = npa_opt(expr, level; goal=:minimise, kw...)
