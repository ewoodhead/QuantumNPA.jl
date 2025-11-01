function sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end

    return matrix
end



"""
    npa_moment(operators)

Construct the NPA moment matrix, given a vector of operators (monomials or
polynomials). The moment matrix returned is a representation of the real part
of the moment matrix, i.e., the moment matrix plus its complex conjugate
divided by two. It is returned in the form of a polynomial with sparse
matrices as the monomials.

# Examples

```julia-repl
julia> gamma = npa_moment([Id, A1, A2, B1, B2]);

julia> gamma[A1]
5×5 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
  ⋅   1.0   ⋅    ⋅    ⋅
 1.0   ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅
  ⋅    ⋅    ⋅    ⋅    ⋅

julia> QuantumNPA.repack(gamma)
5×5 Matrix{Polynomial}:
 Id  A1     A2     B1     B2
 A1  Id     A1 A2  A1 B1  A1 B2
 A2  A1 A2  Id     A2 B1  A2 B2
 B1  A1 B1  A2 B1  Id     B1 B2
 B2  A1 B2  A2 B2  B1 B2  Id
```

```julia-repl
julia> gamma = npa_moment([Id, A1 + A2, A1 - A2]);

julia> gamma[A1]
3×3 SparseArrays.SparseMatrixCSC{Float64, Int64} with 4 stored entries:
  ⋅   1.0  1.0
 1.0   ⋅    ⋅
 1.0   ⋅    ⋅


julia> gamma[A2]
3×3 SparseArrays.SparseMatrixCSC{Float64, Int64} with 4 stored entries:
   ⋅   1.0  -1.0
  1.0   ⋅     ⋅
 -1.0   ⋅     ⋅

julia> QuantumNPA.repack(gamma)
3×3 Matrix{Polynomial}:
 Id       A1 + A2             A1 - A2
 A1 + A2  2.0 Id + 2.0 A1 A2  0
 A1 - A2  0                   2.0 Id - 2.0 A1 A2
```
"""
function npa_moment(operators::Vector)
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

"""
    npa_moment(source, level)

Constructs the NPA moment matrix at the given level of the hierarchy, taking
all degree 1 monomials appearing in `source` as the level 1 operators.
`source` can be a monomial, polynomial, or collection containing monomials,
polynomials, or sub-collections. The level can be a nonnegative integer or a
string, such as `"1 + A B"`.
"""
npa_moment(source, level) = npa_moment(npa_level(source, level))



"""
    npa2sdp(expr, level; eq=[], ge=[])

Generate the NPA relaxation for a given quantum optimisation problem, at the
indicated level of the NPA hierarchy.

This basically uses the `npa_moment` function to generate the moment matrix
and then calls `npa2sdp` again with the moment matrix as the second argument,
which takes the real part of the problem and then eliminates the equality
constraints by substitution.

# Example

```julia-repl
julia> S = A1*(B1 + B2) + A2*(B1 - B2)
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> (objective, moments) = npa2sdp(S, 1, eq=[A2*B1 - Id]);

julia> objective
Id + A1 B1 + A1 B2 - A2 B2

julia> moments[1][Id]
5×5 SparseArrays.SparseMatrixCSC{Float64, Int64} with 7 stored entries:
 1.0   ⋅    ⋅    ⋅    ⋅
  ⋅   1.0   ⋅    ⋅    ⋅
  ⋅    ⋅   1.0  1.0   ⋅
  ⋅    ⋅   1.0  1.0   ⋅
  ⋅    ⋅    ⋅    ⋅   1.0

julia> QuantumNPA.repack(moments[1])
5×5 Matrix{Polynomial}:
 Id  A1     A2     B1     B2
 A1  Id     A1 A2  A1 B1  A1 B2
 A2  A1 A2  Id     Id     A2 B2
 B1  A1 B1  Id     Id     B1 B2
 B2  A1 B2  A2 B2  B1 B2  Id
```
"""
function npa2sdp(expr, level; eq=[], ge=[])
    moment = npa_moment([expr, eq, ge], level)
    return npa2sdp(expr, moment, eq=eq, ge=ge)
end

"""
    npa2sdp(expr, moment; eq=[], ge=[])

Generate the NPA relaxation for a given quantum optimisation problem.

The main purpose of this function is to eliminate equality constraints from
the input problem. This function also takes a representation of the real part
of the input problem, if it isn't already real-valued.

The return value is a tuple containing 1) the expression to optimise, and 2) a
vector of matrices that will be required to be positive semidefinite. The
first element of the vector is derived from the input moment matrix; any
additional ones are derived from any additional inequality constraints
provided.

# Arguments

- `expr`: the operator whose expectation value we want to optimise.
- `moment`: the moment matrix
- `eq`=[] and `ge=[]`: vectors of operators whose expectation values we want
  to set equal to and lower bound by zero, respectively.
"""
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
# Made more specific (Matrix instead of AbstractMatrix) to avoid ambiguity in Julia 1.12
function LinearAlgebra.dot(A::SparseMatrixCSC,
                           B::Symmetric{<:JuMP._MA.AbstractMutable, <:Matrix})
    acc = zero(eltype(B))

    for j in 1:size(A, 2)
        for k in nzrange(A, j)
            add_to_expression!(acc, nonzeros(A)[k], B[rowvals(A)[k], j])
        end
    end

    return acc
end

function LinearAlgebra.dot(A::Symmetric{<:JuMP._MA.AbstractMutable, <:Matrix},
                           B::SparseMatrixCSC)
    return dot(B, A)
end

# Julia 1.12 compatibility: handle plain matrices of JuMP variables
# Need to be more specific than AbstractMatrix to avoid ambiguity with MutableArithmetics
function LinearAlgebra.dot(A::SparseMatrixCSC,
                           B::Matrix{<:JuMP._MA.AbstractMutable})
    acc = zero(eltype(B))

    for j in 1:size(A, 2)
        for k in nzrange(A, j)
            add_to_expression!(acc, nonzeros(A)[k], B[rowvals(A)[k], j])
        end
    end

    return acc
end

function LinearAlgebra.dot(A::Matrix{<:JuMP._MA.AbstractMutable},
                           B::SparseMatrixCSC)
    return dot(B, A)
end



function maximise_p(sense::Symbol)
    if sense in (:maximise, :maximize, :max)
        return true
    elseif sense in (:minimise, :minimize, :min)
        return false
    else
        @error "Unrecognised sense: $sense."
    end
end

function set_objective!(model::Model, maximise::Bool, objective)
    if maximise
        @objective(model, Max, objective)
    else
        @objective(model, Min, objective)
    end
end

function set_objective!(model::Model, sense::Symbol, objective)
    set_objective!(model, maximise_p(sense), objective)
end



function sdp2jump_d(expr, ineqs;
                    sense=:maximise,
                    solver=default_solver,
                    silent=default_silent)
    maximise = maximise_p(sense)
    s = (maximise ? 1 : -1)

    model = jump_model(solver, silent)

    Zs = [@variable(model, [1:m, 1:n], PSD)
          for (m, n) in size_as_pair.(ineqs)]

    Ids = (ineq[Id] for ineq in ineqs)
    objective = (s*sum(dot(m,z)
                     for (m, z) in zip(Ids, Zs))
                 + expr[Id])

    set_objective!(model, !maximise, objective)

    mons = collect(m for m in monomials(expr, ineqs) if !isidentity(m))

    for m in mons
        c = expr[m]
        Fs = (ineq[m] for ineq in ineqs)
        tr_term = sum(dot(F, Z)
                      for (F, Z) in zip(Fs, Zs))

        @constraint(model, tr_term + s*c == 0)
    end

    return model
end



function add_constraint!(model::Model, constraint)
    if typeof(size(constraint)) === Tuple{Int, Int}
        @constraint(model, constraint in PSDCone())
    else
        @constraint(model, constraint .>= 0)
    end
end


function jump_scalar_expr(expr, mons, v)
    result = AffExpr(expr[Id])

    for m in mons
        add_to_expression!(result, expr[m], v[m])
    end

    return result
end

function jump_array_expr(expr, mons, v)
    result = Array{AffExpr}(undef, size(expr)...)

    cId = expr[Id]

    for k in 1:length(result)
        result[k] = AffExpr(cId[k])
    end

    for m in mons
        vm = v[m]

        for (i, j, c) in zip(findnz(expr[m])...)
            add_to_expression!(result[i, j], c, vm)
        end
    end

    return result
end

function jump_expr(expr, mons, v)
    if isscalar(expr)
        return jump_scalar_expr(expr, mons, v)
    else
        return jump_array_expr(expr, mons, v)
    end
end

function sdp2jump(expr, ineqs;
                  sense=:maximise,
                  solver=default_solver,
                  silent=default_silent)
    model = jump_model(solver, silent)

    mons = filter(!isidentity, monomials(ineqs))

    @variable(model, v[mons])

    objective = jump_expr(expr, mons, v)
    set_objective!(model, sense, objective)

    for ineq in ineqs
        constraint = jump_expr(ineq, mons, v)
        add_constraint!(model, constraint)
    end

    return model
end



function npa2jump_d(expr, level_or_moments; eq=[], ge=[], kw...)
    (expr, moments) = npa2sdp(expr, level_or_moments, eq=eq, ge=ge)
    model = sdp2jump_d(expr, moments; kw...)
    return model
end

function npa2jump(expr, level_or_moments; eq=[], ge=[], kw...)
    (expr, moments) = npa2sdp(expr, level_or_moments, eq=eq, ge=ge)
    model = sdp2jump(expr, moments; kw...)
    return model
end



function npa_opt(expr, level_or_moments; eq=[], ge=[], kw...)
    model = npa2jump_d(expr, level_or_moments, eq=eq, ge=ge; kw...)
    optimize!(model)
    return objective_value(model)
end



npa_max(expr, level; kw...) = npa_opt(expr, level; sense=:maximise, kw...)
npa_min(expr, level; kw...) = npa_opt(expr, level; sense=:minimise, kw...)
