# Examples

Maximise CHSH at level 2 of the hierarchy:
```julia-repl
julia> @dichotomic A1 A2 B1 B2;

julia> S = A1*(B1 + B2) + A2*(B1 - B2)
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> npa_max(S, 2)
2.828227712681755
```

Maximise Svetlichny at level $1 + A B + A C + B C$:
```julia-repl
julia> @dichotomic A[1:2] B[1:2] C[1:2];

julia> E(x,y,z) = A[x]*B[y]*C[z]
E (generic function with 1 method)

julia> S = -E(1,1,1) + E(1,1,2) + E(1,2,1) + E(1,2,2) + E(2,1,1) + E(2,1,2) + E(2,2,1) - E(2,2,2)
-A1 B1 C1 + A1 B1 C2 + A1 B2 C1 + A1 B2 C2 + A2 B1 C1 + A2 B1 C2 + A2 B2 C1 - A2 B2 C2

julia> npa_max(S, "1 + A B + A C + B C")
5.656854315137034
```
(note that the spaces e.g. between A and B are necessary in the string, since
party labels go from A to Z then AA to ZZ then AAA to ZZZ...)

Maximise a modified CHSH at level $1 + A B + A^{2} B$:
```julia-repl
julia> npa_max(0.3 * A1 + 0.6 * A1*(B1 + B2) + A2*(B1 - B2), "1 + A B + A^2 B")
2.3584761820283977
```

You can specify both equality and inequality arguments using the `eq` and
`ge` keyword arguments. These should be lists of operators whose expectation
values you want, respectively, to set to and lower bound by zero. For
example, to maximise `<A1>` subject to `<A1*(B1 + B2)> = 1.4` and `<A2*(B1 -
B2)> = 1.4`:
```julia-repl
julia> npa_max(A1, 2, eq=[A1*(B1 + B2) - 1.4*Id, A2*(B1 - B2) - 1.4*Id])
0.19800616634180992
```
Maximise `<A1 + A2>` subject to `<A1 + 2*A2> <= 1 ` and `<2*A1 + A2> <= 1`:
```julia-repl
julia> npa_max(A1 + A2, 1, ge=[Id - A1 - 2*A2, Id - 2*A1 - A2])
0.6666666597867417
```
Maximise `<A1 + A2>` subject to `<A1> = <A2>` and `<A1 + 2*A2> <= 1 `:
```julia-repl
julia> npa_max(A1 + A2, 1, eq=[A1 - A2], ge=[1 - A1 - 2*A2])
0.666642228695571
```

The above examples all use dichotomic variables, but projectors are also
supported. Here we maximise the CH74 form of CHSH:
```julia-repl
julia> PA11, PA12 = projector(1,1,1:2);

julia> PB11, PB12 = projector(2,1,1:2);

julia> npa_max(-PA11 - PB11 + PA11*(PB11 + PB12) + PA12*(PB11 - PB12), 1)
0.20701116471401693

julia> (sqrt(2) - 1)/2
0.20710678118654757
```

Maximise CGLMP with $d=3$ at level $1 + A B$:
```julia-repl
julia> npa_max(cglmp(3), "1 + A B")
2.914945976226541

julia> 1 + sqrt(11/3)
2.914854215512676
```
This uses a function `cglmp()` already defined in `qnpa.jl` to construct the
CGLMP operator.

Maximise the global guessing probability $P_{\text{guess}}(A1,B1|E)$ in the
CHSH setting using full statistics:
```julia
# Create projectors. The keyword argument 'full=true' means that the
# operator corresponding to the highest-numbered output is directly set to
# the identity minus the other ones. For example,
#
#   PA[2,1] = Id - PA[1,1]
#
# and
#
#   PE[4] = Id - PE[1] - PE[2] - PE[3]
#
# This is meant to make working in the Collins-Gisin projection (which the
# NPA code uses) more convenient.
PA = projector(1, 1:2, 1:2, full=true)
PB = projector(2, 1:2, 1:2, full=true)
PE = projector(5, 1:4, 1, full=true)

# CHSH = 2*sqrt(2) * p
p = 0.9

# Expectation value of G is the probability that Eve correctly guesses
# Alice's and Bob's joint outcome.
G = sum(PA[a,1] * PB[b,1] * PE[2*(a-1) + b]
        for a in 1:2 for b in 1:2)

# Ideal CHSH-violating correlations mixed with noise. N.B., the actual
# constraints imposed are on the expectation values of the operators
# in the array.
constraints = [PA[1,1] - 0.5*Id,
               PA[1,2] - 0.5*Id,
               PB[1,1] - 0.5*Id,
               PB[1,2] - 0.5*Id,
               PA[1,1]*PB[1,1] - 0.25*(1 + p/sqrt(2))*Id,
               PA[1,1]*PB[1,2] - 0.25*(1 + p/sqrt(2))*Id,
               PA[1,2]*PB[1,1] - 0.25*(1 + p/sqrt(2))*Id,
               PA[1,2]*PB[1,2] - 0.25*(1 - p/sqrt(2))*Id]

# This returns about 0.7467 for p = 0.9 at level 2 using the default SCS
# solver.
npa_max(G, 2, eq=constraints)
```

QuantumNPA calls the SCS solver by default (since it doesn't require a
license) to solve the NPA relaxation of a quantum optimisation problem, but a
keyword argument `solver` lets you specify a different one. E.g., solve a
problem using Mosek (which you need a license file to use):
```julia-repl
julia> using MosekTools

julia> npa_max(S, 2, solver=Mosek.Optimizer)
2.82842711211242
```
You can also change the default solver if you don't want to specify it every
time, e.g.,
```julia-repl
julia> set_solver!(Mosek.Optimizer)
```

If you want to construct a JuMP model and solve it separately:
```julia-repl
julia> model = npa2jump(S, "1 + A B", solver=SCS.Optimizer)
A JuMP Model
Minimization problem with:
Variables: 45
Objective function type: AffExpr
`AffExpr`-in-`MathOptInterface.EqualTo{Float64}`: 16 constraints
`Vector{VariableRef}`-in-`MathOptInterface.PositiveSemidefiniteConeTriangle`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: SCS

julia> optimize!(model)
------------------------------------------------------------------
               SCS v3.2.0 - Splitting Conic Solver
        (c) Brendan O'Donoghue, Stanford University, 2012
------------------------------------------------------------------
problem:  variables n: 45, constraints m: 61
cones:    z: primal zero / dual free vars: 16
          s: psd vars: 45, ssize: 1
settings: eps_abs: 1.0e-04, eps_rel: 1.0e-04, eps_infeas: 1.0e-07
          alpha: 1.50, scale: 1.00e-01, adaptive_scale: 1
          max_iters: 100000, normalize: 1, rho_x: 1.00e-06
          acceleration_lookback: 10, acceleration_interval: 10
lin-sys:  sparse-direct
          nnz(A): 81, nnz(P): 0
------------------------------------------------------------------
 iter | pri res | dua res |   gap   |   obj   |  scale  | time (s)
------------------------------------------------------------------
     0| 2.22e+01  1.00e+00  2.00e+02 -9.98e+01  1.00e-01  1.40e-04 
    75| 1.51e-05  5.94e-08  5.35e-08  2.83e+00  1.00e-01  3.10e-03 
------------------------------------------------------------------
status:  solved
timings: total: 3.48e-03s = setup: 3.66e-04s + solve: 3.11e-03s
         lin-sys: 2.47e-04s, cones: 2.58e-03s, accel: 1.57e-05s
------------------------------------------------------------------
objective = 2.828427
------------------------------------------------------------------

julia> objective_value(model)
2.8284273129779325
```
If you call `npa2jump()` without the `solver` keyword argument then a solver
isn't assigned, and you will have to assign one to the JuMP model using
JuMP's `set_optimizer()` function. You can also suppress the output of the
solver by either calling `npa2jump()` with the keyword argument `verbose` set
to false or by using JuMP's `set_silent()` function on the returned JuMP
model.
