# QuantumNPA

Code to do NPA in Julia. In development - names of important functions or
even the entire project could change.

Prerequisites:
```julia
using Pkg; Pkg.add(["Combinatorics", "JuMP", "SCS", "BlockDiagonals"])
```

Then to use or try out:
```julia
include("QuantumNPA.jl");
using .QuantumNPA
```
(The dot in the second line isn't a typo.)

If you are working on the code and want to be able to call internal functions
more conveniently you can instead do
```julia
include("qnpa.jl")
```
This will load every function and global variable into the main module.



## Working examples

Maximise CHSH at level 2 of the hierarchy:
```julia
julia> @dichotomic A1 A2 B1 B2;

julia> S = A1*(B1 + B2) + A2*(B1 - B2)
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> npa_max(S, 2)
2.828227712681755
```

Maximise Svetlichny at level 1 + A B + A C + B C:
```julia
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

Maximise a modified CHSH at level 1 + A B + A^2 B:
```julia
julia> npa_max(0.3 * A1 + 0.6 * A1*(B1 + B2) + A2*(B1 - B2), "1 + A B + A^2 B")
2.3584761820283977
```

You can specify both equality and inequality arguments using the `eq` and
`ge` keyword arguments. These should be lists of operators whose expectation
values you want, respectively, to set to and lower bound by zero. For
example, to maximise `<A1>` subject to `<A1*(B1 + B2)> = 1.4` and `<A2*(B1 -
B2)> = 1.4`:
```julia
julia> npa_max(A1, 2, eq=[A1*(B1 + B2) - 1.4, A2*(B1 - B2) - 1.4])
0.19800616634180992
```
Maximise `<A1 + A2>` subject to `<A1 + 2*A2> <= 1 ` and `<2*A1 + A2> <= 1`:
```julia
julia> npa_max(A1 + A2, 1, ge=[1 - A1 - 2*A2, 1 - 2*A1 - A2])
0.6666666597867417
```
Maximise `<A1 + A2>` subject to `<A1> = <A2>` and `<A1 + 2*A2> <= 1 `:
```julia
julia> npa_max(A1 + A2, 1, eq=[A1 - A2], ge=[1 - A1 - 2*A2])
0.666642228695571
```

The above examples all use dichotomic variables, but projectors are also
supported. Here we maximise the CH74 form of CHSH:
```julia
julia> PA11, PA12 = projector(1,1,1:2);

julia> PB11, PB12 = projector(2,1,1:2);

julia> npa_max(-PA11 - PB11 + PA11*(PB11 + PB12) + PA12*(PB11 - PB12), 1)
0.20701116471401693

julia> (sqrt(2) - 1)/2
0.20710678118654757
```

Maximise CGLMP with d=3 at level 1 + A B:
```julia
julia> npa_max(cglmp(3), "1 + A B")
2.914945976226541

julia> 1 + sqrt(11/3)
2.914854215512676
```
This uses a function `cglmp()` already defined in `qnpa.jl` to construct the
CGLMP operator.

Maximise the global guessing probability Pg(A1,B1|E) in the CHSH setting
using full statistics:
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
constraints = [PA[1,1] - 0.5,
               PA[1,2] - 0.5,
               PB[1,1] - 0.5,
               PB[1,2] - 0.5,
               PA[1,1]*PB[1,1] - 0.25*(1 + p/sqrt(2)),
               PA[1,1]*PB[1,2] - 0.25*(1 + p/sqrt(2)),
               PA[1,2]*PB[1,1] - 0.25*(1 + p/sqrt(2)),
               PA[1,2]*PB[1,2] - 0.25*(1 - p/sqrt(2))]

# This returns about 0.7467 for p = 0.9 at level 2 using the default SCS
# solver.
npa_max(G, 2, eq=constraints)
```

QuantumNPA calls the SCS solver by default (since it doesn't require a
license) to solve the NPA relaxation of a quantum optimisation problem, but a
keyword argument `solver` lets you specify a different one. E.g., solve a
problem using Mosek (which you need a license file to use):
```julia
julia> using MosekTools

julia> npa_max(S, 2, solver=Mosek.Optimizer)
2.82842711211242
```
You can also change the default solver if you don't want to specify it every
time, e.g.,
```julia
julia> set_solver!(Mosek.Optimizer)
```

If you want to construct a JuMP model and solve it separately:
```julia
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



## Basic features

We can do arithmetic with and take conjugates of some different types of
operators that we associate to different parties. At the moment:
- dichotomic,
- fourier,
- projector,
- unitary,
- zbff (operators for Brown-Fawzi-Fawzi method).

The identity is represented by a variable `Id` that is predefined.
```julia
julia> Id
Id

julia> projector(1, 2, 3)
PA2|3

julia> PA = projector(1, 1:2, 1:2);

julia> PB = projector(2, 1:2, 1:2);

julia> PA[1,1]
PA1|1

julia> PA[1,1]*PB[1,1]
PA1|1 PB1|1

julia> PB[1,1]*PA[1,1]
PA1|1 PB1|1

julia> PA[1,1]*PA[2,1]
0

julia> PA[1,1]*PA[1,1]
PA1|1

julia> UA = unitary(1, 1:3);

julia> V = UA[1]*conj(UA[2])*UA[3]
UA1 UA*2 UA3

julia> P = PA[1,1] + V
PA1|1 + UA1 UA*2 UA3

julia> conj(P)
PA1|1 + UA*3 UA2 UA*1

julia> P*P
PA1|1 + PA1|1 UA1 UA*2 UA3 + UA1 UA*2 UA3 PA1|1 + UA1 UA*2 UA3 UA1 UA*2 UA3

julia> conj(P)*P
Id + PA1|1 + PA1|1 UA1 UA*2 UA3 + UA*3 UA2 UA*1 PA1|1

julia> Q = Id + V*PA[1,1]
Id + UA1 UA*2 UA3 PA1|1

julia> Q*Q
Id + 2 UA1 UA*2 UA3 PA1|1 + UA1 UA*2 UA3 PA1|1 UA1 UA*2 UA3 PA1|1

julia> conj(Q)*Q
Id + PA1|1 + PA1|1 UA*3 UA2 UA*1 + UA1 UA*2 UA3 PA1|1

julia> ZE = zbff(5, 1:2);

julia> R = PA[1,1]*PB[2,2]*ZE[1]*ZE[2]
PA1|1 PB2|2 ZE1 ZE2

julia> conj(R)
PA1|1 PB2|2 ZE*2 ZE*1

julia> FA = fourier(1, 9, 1, 5)
A9^1

julia> FA^0
Id

julia> FA^3
A9^3

julia> conj(FA^3)
A9^2

julia> FA^5
Id

julia> FA^6
A9^1

julia> conj(FA^3)*FA^3
Id

julia> FA*FA
A9^2

julia> conj((FA*FA)^4)
A9^2

julia> A1, A2 = dichotomic(1, 1:2);

julia> B1, B2 = dichotomic(2, 1:2);

julia> S = A1*(B1 + B2) + A2*(B1 - B2)
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> S^2
4 Id - A1 A2 B1 B2 + A1 A2 B2 B1 + A2 A1 B1 B2 - A2 A1 B2 B1
```

The functions that create monomials and their parameters are
```julia
dichotomic(party, input)
fourier(party, input, power, d)
projector(party, output, input, full=false)
unitary(party, index, conj=false)
zbff(party, index, conj=false)
```
In these:
- `party` is either:
  - a strictly positive integer, e.g., `1`,
  - a vector of strictly positive integers in strictly increasing order,
    e.g., `[1, 3, 4]`, representing an operator associated to more than one
    party,
  - an uppercase alphabetic character, e.g. `'A'`, or
  - a string of alphabetic characters corresponding to parties in increasing
    order separated by underscores, e.g., `"A"` (same as party number `1`),
    `"AB"` (party `28`), `"A_B"` (party vector `[1, 2]`), `"A_BC"` (party
    vector `[1, 55]`).
  Parties are always converted to and stored internally in monomials in the
  vector-of-integers form. Operators associated to different parties are
  considered to commute if and only if the intersection of the party vectors
  is empty.
- The parameters called `input`, `output`, and `index` can be either integers
  or arrays or ranges of integers.
- The parameter `conj` is optional and defaults to `false` if it is omitted.
- For projectors, if you give a range of inputs you can also give a value for
  a fourth parameter `full`, which defaults to `false`. Setting it to `true`
  indicates that you intend for the range of outputs to represent the full
  set of measurement outcomes. In that case, in place of the last projector
  you are given the identity minus the sum of all the preceding projectors.

A few examples illustrating different ways of calling the `projector()`
function:
```julia
julia> projector(1, 1:2, 1:2)
2×2 Array{Monomial,2}:
 PA1|1  PA1|2
 PA2|1  PA2|2

julia> projector(1, 1:3, 1:2, full=true)
3×2 Array{Any,2}:
 PA1|1               PA1|2
 PA2|1               PA2|2
 Id - PA1|1 - PA2|1  Id - PA1|2 - PA2|2

julia> zbff(1, 1:3)
3-element Array{Monomial,1}:
 ZA1
 ZA2
 ZA3
```

Examples illustrating commutation relations with dichotomic operators:
```julia
julia> A1, A2 = dichotomic(1, 1:2);

julia> B1, B2 = dichotomic(2, 1:2);

julia> A_C1 = dichotomic([1,3], 1);

julia> A1*B1
A1 B1

julia> B1*A1
A1 B1

julia> A1*A_C1
A1 A_C1

julia> A_C1*A1
A_C1 A1

julia> A1*B1*A_C1
A1 A_C1 B1

julia> A_C1*A1*B1
A_C1 A1 B1

julia> A_C1*B1*A1
A_C1 A1 B1

julia> B1*A_C1*A1
A_C1 A1 B1

julia> A1*B1*A_C1*A2*B2
A1 A_C1 A2 B1 B2

julia> A1*B1*A_C1*A1*B1
A1 A_C1 A1
```

I am working on writing macros to automatically create variables using
"standard" names. At the moment you can do, e.g., this to create some
dichotomic variables:
```julia
@dichotomic A1 A2 B1 B2 C[1:3]
```
The above macro invocation does essentially the same as running the
following:
```julia
A1 = dichotomic(1, 1)
A2 = dichotomic(1, 2)
B1 = dichotomic(2, 1)
B2 = dichotomic(2, 2)
C = dichotomic(3, 1:3)
```

There are no special relations (at least, at the moment) between the
different types of operators, so you shouldn't, for example, mix projectors
and dichotomic operators unless you consider them to be unrelated to each
other:
```julia
julia> dichotomic(1, 1) * projector(1, 1, 1)
A1 PA1|1

julia> dichotomic(1, 1) - (2*projector(1, 1, 1) - Id)
Id + A1 - 2 PA1|1
```

Finally, note that operators are objects that can be manipulated in the same
sorts of ways as other types of objects in Julia, such as putting them in
arrays or other data structures. For example, `dichotomic(p, 1:n)` returns a
one-dimensional array of dichotomic operators, which we can then use in
vector expressions such as:
```julia
julia> A = dichotomic(1, 1:2)
2-element Vector{Monomial}:
 A1
 A2

julia> B = dichotomic(2, 1:2)
2-element Vector{Monomial}:
 B1
 B2

julia> M = [1 1; 1 -1]
2×2 Matrix{Int64}:
 1   1
 1  -1

julia> A'*M*B
A1 B1 + A1 B2 + A2 B1 - A2 B2
```



## Analysing, modifying, and deconstructing operators

Monomials and polynomials are objects of different types, although a
polynomial consisting of a single monomial multiplied by 1 is printed the
same as a monomial:
```julia
julia> P = projector(1, 1, 1)
PA1|1

julia> typeof(P)
Monomial

julia> Q = 1*P
PA1|1

julia> typeof(Q)
Polynomial

julia> S
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> typeof(S)
Polynomial
```
If you need to ensure a given object is a polynomial you can "promote" it by
calling `Polynomial()` on it. This does nothing if the argument is already a
polynomial:
```julia
julia> x = Polynomial(1)
Id

julia> y = Polynomial(Id)
Id

julia> z = Polynomial(S)
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> typeof.([x, y, z])
3-element Array{DataType,1}:
 Polynomial
 Polynomial
 Polynomial
```
Note that, in the last case, the polynomial returned is the same as (and not
a copy of) the original, which means that modifying `z` here will modify `S`
since they are the same object:
```julia
julia> z === S
true

julia> z[A1] = 7
7

julia> S
7 A1 + A1 B1 + A1 B2 + A2 B1 - A2 B2
```
If you want to create a copy of a polynomial that you can safely modify
without changing the original you can call the `copy()` function to do this:
```julia
julia> S = A1*(B1 + B2) + A2*(B1 - B2)
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> z = copy(S)
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> z === S
false

julia> z[A1] = 7
7

julia> z
7 A1 + A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> S
A1 B1 + A1 B2 + A2 B1 - A2 B2
```

As the two above examples suggest, you can access and/or modify the
coefficient associated to a given monomial using `[]`:
```julia
julia> S[A1*B1]
1

julia> S[A1]
0
```
You can get all the monomials in a polynomial by calling the `monomials()`
function on it:
```julia
julia> monomials(S)
Base.KeySet for a Dict{Monomial,Number} with 4 entries. Keys:
  A1 B2
  A1 B1
  A2 B1
  A2 B2
```
Polynomials will also act as iterators over pairs of their nonzero
coefficients and monomials in contexts where an iterator is expected:
```julia
julia> collect(S)
4-element Vector{Any}:
 Pair{Number,Monomial}(-1, A2 B2)
 Pair{Number,Monomial}(1, A1 B2)
 Pair{Number,Monomial}(1, A2 B1)
 Pair{Number,Monomial}(1, A1 B1)

julia> for (c, m) in S
           @printf "%s  =>  %2d\n" m c
       end
A2 B1  =>   1
A1 B2  =>   1
A2 B2  =>  -1
A1 B1  =>   1
```
If you want to iterate over the monomials in lexicographical order you can
just call `sort()` on the polynomial first:
```julia
julia> for (c, m) in sort(S)
           @printf "%s  =>  %2d\n" m c
       end
A1 B1  =>   1
A1 B2  =>   1
A2 B1  =>   1
A2 B2  =>  -1
```

In order to help analyse a problem, there is a function `operators()` that
can find and return all the individual (order 1) operators in one or more
monomials and polynomials or collections of such operators. For instance, if
we wanted to maximise the local guessing probability in the CHSH setting
using full statistics we might represent the problem by an objective
polynomial and a list of constraint polynomials whose expectation values we
want to set to zero:
```julia
A1, A2 = dichotomic(1, 1:2)
B1, B2 = dichotomic(2, 1:2)
E1 = dichotomic(5, 1)

objective = (1 + A1*E1)/2

constraints = [A1, A2, B1, B2,
               A1*B1 - 0.7,
               A1*B2 - 0.7,
               A2*B1 - 0.7,
               A2*B2 + 0.7]
```
Assuming these variable definitions, we can use the `operators()` function to
immediately find all the level-one operators in the problem:
```julia
julia> operators(objective, constraints)
Set{Monomial} with 5 elements:
  A2
  A1
  B1
  B2
  E1
```
`operators()` can optionally take a keyword argument `by_party` which is set
to `false` by default. Setting it to `true` groups the level-one operators by
party and returns a dictionary of the parties and operators associated to
those parties:
```julia
julia> operators(objective, constraints, by_party=true)
Dict{Integer,Set{Monomial}} with 3 entries:
  2 => Set(Monomial[B1, B2])
  5 => Set(Monomial[E1])
  1 => Set(Monomial[A2, A1])
```
This is useful for constructing operators at intermediate levels of the
NPA hierarchy like "1 + A B + A E + B E".



## NPA example

This short example finds what operators appear and where in the NPA moment
matrix at level 2 for the CHSH problem. It covers only the upper triangular
part and treats monomials and their conjugates as the same.
```julia
A1, A2 = dichotomic(1, 1:2)
B1, B2 = dichotomic(2, 1:2)
ops1 = [Id, A1, A2, B1, B2]
ops2 = sort(Set(O1*O2 for O1 in ops1 for O2 in ops1))

indices = Dict()

indexed_ops = collect(enumerate(ops2))

for (i, x) in indexed_ops
    for (j, y) in indexed_ops[i:end]
        m = conj(x)*y
        m = min(m, conj(m))

        if haskey(indices, m)
            push!(indices[m], (i, j))
        else
            indices[m] = [(i, j)]
        end
    end
end
```
This gives:
```julia
julia> for (m, l) in sort!(collect(indices), by=first)
           @printf "%11s  =>  %s\n" m l
       end
         Id  =>  [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (9, 9), (10, 10), (11, 11), (12, 12), (13, 13)]
         A1  =>  [(1, 2), (3, 7), (4, 8), (5, 9)]
         A2  =>  [(1, 3), (2, 6), (4, 10), (5, 11)]
         B1  =>  [(1, 4), (2, 8), (3, 10), (5, 13)]
         B2  =>  [(1, 5), (2, 9), (3, 11), (4, 12)]
      A1 A2  =>  [(1, 6), (1, 7), (2, 3), (8, 10), (9, 11)]
      A1 B1  =>  [(1, 8), (2, 4), (7, 10), (9, 13)]
      A1 B2  =>  [(1, 9), (2, 5), (7, 11), (8, 12)]
      A2 B1  =>  [(1, 10), (3, 4), (6, 8), (11, 13)]
      A2 B2  =>  [(1, 11), (3, 5), (6, 9), (10, 12)]
      B1 B2  =>  [(1, 12), (1, 13), (4, 5), (8, 9), (10, 11)]
   A1 A2 A1  =>  [(2, 7)]
   A2 A1 A2  =>  [(3, 6)]
   A1 A2 B1  =>  [(2, 10), (3, 8), (4, 6), (4, 7)]
   A1 A2 B2  =>  [(2, 11), (3, 9), (5, 6), (5, 7)]
   A1 B1 B2  =>  [(2, 12), (2, 13), (4, 9), (5, 8)]
   A2 B1 B2  =>  [(3, 12), (3, 13), (4, 11), (5, 10)]
   B1 B2 B1  =>  [(4, 13)]
   B2 B1 B2  =>  [(5, 12)]
A1 A2 A1 A2  =>  [(6, 7)]
A1 A2 A1 B1  =>  [(7, 8)]
A1 A2 A1 B2  =>  [(7, 9)]
A2 A1 A2 B1  =>  [(6, 10)]
A2 A1 A2 B2  =>  [(6, 11)]
A1 A2 B1 B2  =>  [(6, 13), (7, 12), (8, 11)]
A1 A2 B2 B1  =>  [(6, 12), (7, 13), (9, 10)]
A1 B1 B2 B1  =>  [(8, 13)]
A1 B2 B1 B2  =>  [(9, 12)]
A2 B1 B2 B1  =>  [(10, 13)]
A2 B2 B1 B2  =>  [(11, 12)]
B1 B2 B1 B2  =>  [(12, 13)]
```

The example above uses `min(m, conj(m))` to find which of `m` or its
conjugate comes first lexicographically. It works because comparisons between
monomials are defined:
```julia
julia> A1 == A1
true

julia> A1 < A1
false

julia> A1 < A2
true

julia> A2 < A1*A2
true
```
`sort` used above works for the same reason. `==` and `!=` (but not the
inequalities) can also be used to compare monomials with polynomials or
polynomials with each other:
```julia
julia> A1 == 1*A1
true

julia> A1 < 1*A1
ERROR: MethodError: no method matching isless(::Monomial, ::Polynomial)
[...]
```



## Comments on implementation

This last section gives some details about how operators are
implemented. This is mainly of interest to people who want to better
understand what the code does and how it works or want to modify it.

The basic idea behind this whole library is that the NPA hierarchy method
itself is pretty straightforward if the types of operators you need are
supported. One then just has to compute all the operator products appearing
in the upper triangular part of the moment matrix and check for relations
between the results that translate to constraints on the moment matrix. So
the main thing this library aims to do is add support to Julia for doing
arithmetic and automatically simplifying products of certain types of
operators commonly encountered in quantum optimisation problems that can be
handled with the NPA hierarchy.

There are three main types of operator defined in the code. They are:
- The abstract `Operator` type, defined in `src/ops_primitive.jl`, which has
  several concrete subtypes (`Dichotomic`, etc.) defined in
  `src/ops_predefined.jl` corresponding to the different types of supported
  operator. These types and objects of these types are used internally. They
  are not meant to be created or manipulated directly by the user.
  
- `Monomial`, defined in `src/ops_monomial.jl`, to represent products of
  operators associated to different parties.

- `Polynomial`, defined in `src/ops_polynomial.jl`, to represent linear
  combinations of monomials.


### `Operator`

`Operator` is an abstract type. This essentially means it is a collective
name for several concrete types that, throughout the code, are assumed to
have certain common properties making them interchangeable to a significant
extent. In particular, they can be passed as arguments to certain functions
such as `conj()` and `Base.:*`. The concrete subtypes, such as `Dichotomic`,
are structs representing the different types of operators supported. What
fields they have depends on the type. For example, `Dichotomic` objects have
one field `index` for the index, `Projector`s have two for the output and
input, and unitaries have an integer `index` and boolean `conj` field
tracking whether they are conjugated or not.

Objects of the different basic operator types are structs containing data
about them. They are created by calling their constructors, which are
functions with the same (capitalised) names as the types. For example:
```julia
julia> x = Dichotomic(2)
/2

julia> fieldnames(Dichotomic)
(:input,)

julia> x.input
2

julia> y = Projector(2, 3)
P2|3

julia> fieldnames(Projector)
(:output, :input)

julia> (y.output, y.input)
(2, 3)

julia> u = Unitary(3, true)
U*3

julia> fieldnames(Unitary)
(:index, :conj)

julia> (u.index, u.conj)
(3, true)
```
(Note that these constructors are not exported by the `QuantumNPA` module;
you either need to prefix their names, e.g., `QuantumNPA.Dichotomic`, or load
the codebase with `include("qnpa.jl")` to call them.)

The file `src/ops_primitive.jl` defines default implementations of certain
key functions that have to be supported for all operators, which can then be
overriden when the default behaviour is not correct. For example, `conj()` is
defined generically for operators to just return the original operator:
```julia
Base.conj(o::Operator) = o
```
This is fine for operators that are meant to be Hermitian:
```julia
julia> cx = conj(x)
/2

julia> cx === x
true
```
But non-Hermitian types have conjugates that have to be determined in
different ways depending on the type and therefore need their own specialised
versions of `conj()`. For instance, `Unitary` objects have a `conj` field and
a version of `conj()` specialised for unitaries just toggles it. Its
definition is generated by a macro and amounts to
```julia
Base.conj(u::Unitary) = Unitary(u.index, !u.conj)
```
It is this version that gets called if `conj()` is called with a unitary
argument:
```julia
julia> u.conj
true

julia> v = conj(u)
U3

julia> v.conj
false
```

How to multiply operators and, especially, when the product of two operators
can be simplified is determined by a generic and specialised versions of
`Base.:*`. The generic version is:
```julia
Base.:*(x::Operator, y::Operator) = (1, [x, y])
```
This just means that the default rule is to concatenate operators that are
multiplied together, represented by the list `[x, y]`. The number `1` above
is a multiplicative coefficient. This allows for the possibility of scaling
factors appearing in multiplications in the future (such as something like an
unnormalised projector that squares to a multiple of itself). Currently there
are no such operators defined in the codebase and the coefficient is always
zero or one.

The implementation above is not sufficient for most of the operator types
and, in practice, the generic multiplication rule is usually fallen back on
to multiply operators of different types, e.g.,
```julia
julia> x*y
(1, Operator[/2, P2|3])
```
while multiplication of operators of the same type is usually handled by
specialised versions. Among these, `Dichotomic` objects have among the
simplest nontrivial multiplication rules: the product of two dichotomic
operators is the identity if their `input` fields are the same, otherwise
they just concatenate. The implementation is:
```julia
function Base.:*(x::Dichotomic, y::Dichotomic)
    return (x.input == y.input) ? (1, []) : (1, [x, y])
end
```
with the empty vector `[]` used to represent the identity.

Products are represented by vectors (one-dimensional arrays) of operators. A
function `join_ops()`, defined near the beginning of `src/ops_primitive.jl`,
determines how to multiply products. It takes two vectors of operators,
`opsx` and `opsy`, and basically takes out and multiplies the last element of
`opsx` and `opsy`, repeats this if the result is `[]` (representing the
identity), and then returns the concatenation of the remaining elements of
`opsx`, the last non-empty product, and the remaining `opsy`. An example
invocation:
```julia
julia> p = [u, x]
2-element Vector{Operator}:
 U*3
 /2

julia> q = [x, y]
2-element Vector{Operator}:
 /2
 P2|3

julia> join_ops(p, q)
(1, Operator[U*3, P2|3])
```



### `Monomial`

`Monomial`s are the most basic type of operator that are meant to be created
and manipulated in normal use of QuantumNPA. They represent products of
operators grouped into different parties. They are structs whose only field,
`word`, contains a vector `[(p1, ops1), (p2, ops2), ...]` of pairs of party
vectors `p1`, `p2`, etc. and vectors of operators `ops1`, `ops2` (of the type
used by the `join_ops()` function described above) associated to those
parties. Thus, we can look at the contents of a `Monomial` by accessing its
`word` field:
```julia
julia> (A1, A2) = dichotomic(1, 1:2);

julia> PB11 = projector(2, 1, 1);

julia> (UE1, UE2) = unitary(5, 1:2);

julia> M = A1*A2*PB11*UE1*UE2
A1 A2 PB1|1 UE1 UE2

julia> M.word
3-element Vector{Tuple{Vector{Int64}, Vector{Operator}}}:
 ([1], [/1, /2])
 ([2], [P1|1])
 ([5], [U1, U2])
```

Party vectors (on the left of the example above) are vectors of integers
representing which party or parties a group of operators are associated
to. Valid party vectors are vectors of integers, such as `[1, 2, 4]`, in
which all the integers are in strictly increasing order and the first (and
smallest) integer is at least one. One party vector `p` is considered to
lexicographically precede another `q` if `p < q` returns `true`. Often, they
will just contain a single party, e.g. , `[1]`, but this isn't
required. Operators associated to different parties are taken to commute if
the intersection of the party vectors is empty. Thus `[1]` commutes with
`[2]` and `[1,3]` commutes with `[2,4]`, but `[1,2]` does not commute with,
for example, `[2]` or `[2,3]`

`Monomial` objects are meant to represent monomials in a certain reduced
canonical form. A monomial is considered in correctly reduced form if:

1. The party vectors are valid and appear in lexicographic order as much as
   commutation relations between them allow. This basically means that if a
   party vector `p` is immediately followed by a party vector `q` then at
   least one of `p < q` and `intersect(p, q) != []` should be `true`.

2. The vectors of operators are nonempty and reduced as much as possible. For
   example, a valid vector should not contain the same dichotomic operator
   twice, or a unitary and its conjucate, directly following one another.

A few key functions, particularly `Base.:*()`, `conj()`, and `adjoint()`, are
responsible for maintaining these conventions, i.e., they should return
monomials in the above-described canonical form assuming their inputs are in
canonical form. So the recommended way to build monomials in most cases is to
start with monomials containing a single operator and multiply them to
construct longer monomials. A monomial containing just one operator can be
created by calling the `Monomial` function with a party vector and operator
as arguments, e.g.,
```julia
julia> M = Monomial([1], Dichotomic(2))
A2

julia> M.word
1-element Vector{Tuple{Vector{Int64}, Vector{Operator}}}:
 ([1], [/2])
```
It is also possible to construct a `Monomial` by calling `Monomial()` with an
array of pairs of party vectors and vectors of operators as an argument, but
this way isn't very readable and makes it easy to generate invalid monomials,
and so should be avoided.



### `Polynomial`

`Polynomial` objects represent linear combinations of monomials, such as
`A1 + 2 A2 B1`. They have a single field, `terms`, which is a dictionary
mapping monomials to coefficients:
```julia
julia> P = 3*Id + 2*A1*A2
3 Id + 2 A1 A2

julia> P.terms
Dict{Monomial, Number} with 2 entries:
  Id    => 3
  A1 A2 => 2
```
Only terms with nonzero coeffients should be stored. The arithmetic functions
and indexed assignment (`setindex!()`, which makes `p[p] = c` work) don't
create or remove pairs for which the coefficient is zero. So setting a term
to zero deletes it from the dictionary:
```julia
julia> P
3 Id + 2 A1 A2

julia> P[A1*A2] = 0
0

julia> P
3 Id

julia> P.terms
Dict{Monomial, Number} with 1 entry:
  Id => 3
```
The zero polynomial is represented by an empty dictionary.

There are a few different versions of the `Polynomial()` constructor. The
basic one takes a dictionary mapping monomials to coefficients and simply
uses that as the `terms` field. This is only meant to be used internally, and
with care, since it can be used to create "invalid" polynomials that break
assumptions made elsewhere in the library:
```julia
julia> Q = Polynomial(Dict(A1*A2 => 0))
0 A1 A2

julia> Q == 0
false
```
Other versions create a polynomial with no input argument (returns the zero
polynomial), or a number, a monomial, a number and monomial, or a polynomial:
```julia
julia> Polynomial()
0

julia> Polynomial(3)
3 Id

julia> Polynomial(A1)
A1

julia> Polynomial(3, A1)
3 A1

julia> Polynomial(P)
3 Id
```
As mentioned above, the latter just returns the polynomial given as input.
The fourth one is used to define multiplication of a number by a monomial in
`src/ops_polynomial.jl`:
```julia
Base.:*(x::Number, y::Monomial) = Polynomial(x, y)
Base.:*(x::Monomial, y::Number) = Polynomial(y, x)
```
