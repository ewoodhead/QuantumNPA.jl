# bff-npa

Code to do NPA needed to use the Brown-Fawzi-Fawzi method.

At the moment: `bff.jl` contains code in development, `ops.jl` is the older
code that only does NPA with projectors and shouldn't be used (except to
copy/adapt some of the code in it and move to `bff.jl`).

Use `bff.jl` like this:
```julia
using Printf
using Base.Iterators
include("bff.jl")
```


## Basic features

At the moment it can do arithmetic with projectors and those Z operators
in the Brown-Fawzi-Fawzi paper, e.g.,
```julia
julia> Id
 Id

julia> projector(1,1,1)
 PA1|1

julia> conj(projector(1,1,1))
 PA1|1

julia> projector(1,1,1)*projector(2,1,1)
 PA1|1 PB1|1

julia> projector(1,1,1)*projector(1,2,1)
0

julia> projector(1,1,1)*projector(1,2,2)
 PA1|1 PA2|2

julia> conj(projector(1,1,1)*projector(1,2,2))
 PA2|2 PA1|1

julia> zbff(1,3)
 ZA3

julia> zbff(1,3,true)
 ZA*3

julia> conj(zbff(1,3))
 ZA*3

julia> conj(zbff(1,3,true))
 ZA3

julia> p = projector(1,1,1)*projector(2,3,4)*zbff(5,1)*zbff(5,2)
 PA1|1 PB3|4 ZE1 ZE2

julia> conj(p)
 PA1|1 PB3|4 ZE*2 ZE*1

julia> p*p
 PA1|1 PB3|4 ZE1 ZE2 ZE1 ZE2

julia> conj(p)*p
 PA1|1 PB3|4 ZE*2 ZE*1 ZE1 ZE2

julia> p*conj(p)
 PA1|1 PB3|4 ZE1 ZE2 ZE*2 ZE*1

julia> diop(1,1)
 + (-1) Id + (2) PA1|1
```

Syntax is
```julia
projector(party, output, input)
zbff(party, index)
zbff(party, index, conj)
diop(party, input)
```
Party numbers start from 1. `diop(p, i)` makes a dichotomic operator; it does
the same thing as `2*projector(p, 1, i) - Id`. The parameters `output` and
`input` (for projectors and diops) and `index` for the Z operators can also
be ranges or lists of integers. The functions return arrays of operators in
these cases:
```julia
julia> projector(1,1:2,1:2)
2×2 Array{Monomial,2}:
  PA1|1   PA1|2
  PA2|1   PA2|2

julia> zbff(1,1:3)
3-element Array{Monomial,1}:
  ZA1
  ZA2
  ZA3

julia> diop(1,1:3)
3-element Array{Polynomial,1}:
  + (-1) Id + (2) PA1|1
  + (-1) Id + (2) PA1|2
  + (-1) Id + (2) PA1|3
```

It is also possible to do general arithmetic with operators, for example:
```julia
julia> A = diop(1,1:2); B = diop(2,1:2);

julia> S = A[1]*(B[1] + B[2]) + A[2]*(B[1] - B[2])
 + (2) Id + (-4) PA1|1 + (-4) PB1|1 + (4) PA1|1 PB1|1 + (4) PA1|1 PB1|2 + (4) PA1|2 PB1|1 + (-4) PA1|2 PB1|2
```
Note that monomials and polynomials are different types, and it is possible
to loop over the monomials and (nonzero) coefficients in a polynomial:
```julia
julia> projector(1,1,1)
 PA1|1

julia> typeof(projector(1,1,1))
Monomial

julia> 1*projector(1,1,1)
 + (1) PA1|1

julia> typeof(1*projector(1,1,1))
Polynomial

julia> A[1]
 + (-1) Id + (2) PA1|1

julia> typeof(A[1])
Polynomial

julia> typeof(S)
Polynomial

julia> for (m, c) in S
           @printf "%12s  =>  %2d\n" m c
       end
 PA1|1 PB1|2  =>   4
       PA1|1  =>  -4
 PA1|2 PB1|1  =>   4
 PA1|2 PB1|2  =>  -4
 PA1|1 PB1|1  =>   4
       PB1|1  =>  -4
          Id  =>   2

julia> for (m, c) in sort(S)
           @printf "%12s  =>  %2d\n" m c
       end
          Id  =>   2
       PA1|1  =>  -4
       PB1|1  =>  -4
 PA1|1 PB1|1  =>   4
 PA1|1 PB1|2  =>   4
 PA1|2 PB1|1  =>   4
 PA1|2 PB1|2  =>  -4
```


## NPA example

This short example finds what operators appear and where in the NPA moment
matrix at level 2 for the CHSH problem. It covers only the upper triangular
part and treats monomials and their conjugates as the same.
```julia
PA, PB = projector(1,1,1:2), projector(2,1,1:2)
ops1 = [PA[1], PA[2], PB[1], PB[2]]
ops2 = sort([pA*pB for pA in PA for pB in PB])
ops = vcat([Id], ops1, ops2)

indices = Dict()

for (i, x) in enumerate(ops)
    for j in i:length(ops)
        y = ops[j]
        m = conj(x)*y
        m = min(m, conj(m))

        if m == 0
            continue
        elseif haskey(indices, m)
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
           @printf "%24s  =>  %s\n" m l
       end
                      Id  =>  [(1, 1)]
                   PA1|1  =>  [(1, 2), (2, 2)]
                   PA1|2  =>  [(1, 3), (3, 3)]
                   PB1|1  =>  [(1, 4), (4, 4)]
                   PB1|2  =>  [(1, 5), (5, 5)]
             PA1|1 PA1|2  =>  [(2, 3)]
             PA1|1 PB1|1  =>  [(1, 6), (2, 4), (2, 6), (4, 6), (6, 6)]
             PA1|1 PB1|2  =>  [(1, 7), (2, 5), (2, 7), (5, 7), (7, 7)]
             PA1|2 PB1|1  =>  [(1, 8), (3, 4), (3, 8), (4, 8), (8, 8)]
             PA1|2 PB1|2  =>  [(1, 9), (3, 5), (3, 9), (5, 9), (9, 9)]
             PB1|1 PB1|2  =>  [(4, 5)]
       PA1|1 PA1|2 PB1|1  =>  [(2, 8), (3, 6), (6, 8)]
       PA1|1 PA1|2 PB1|2  =>  [(2, 9), (3, 7), (7, 9)]
       PA1|1 PB1|1 PB1|2  =>  [(4, 7), (5, 6), (6, 7)]
       PA1|2 PB1|1 PB1|2  =>  [(4, 9), (5, 8), (8, 9)]
 PA1|1 PA1|2 PB1|1 PB1|2  =>  [(6, 9)]
 PA1|1 PA1|2 PB1|2 PB1|1  =>  [(7, 8)]
```

The example above uses `min(m, conj(m))` to find which of `m` or its
conjugate comes first lexicographically. It works because comparisons between
monomials are defined:
```julia
julia> PA[1] == PA[1]*PA[1]
true

julia> PA[1] < PA[1]
false

julia> PA[1] < PA[2]
true
```
`sort` used above works for the same reason. `==` and `!=` (but not the
inequalities) can also be used to compare polynomials.


## Internal details

The way a list of operators are joined to multiply them is determined at the
moment by a function `join_ops` near the beginning of the file `bff.jl`. This
is what would need to be generalised if we wanted to support multiplication
of slightly more general types of operators. For example, if we introduced a
unitary type with `U* U = 1`, the function `join_ops` as it is at the moment
would simplify `U* U* U U` to `U* U` but not all the way to `1`.

Associating operators in groups to parties is handled by the `Monomial`
type. At the moment it just contains a list
```julia
word = [(p1, ops1), (p2, ops2), ...]
```
of parties `p1`, `p2`, etc. and lists of operators `ops1`, `ops2`,
etc. associated to those parties. For example,
```julia
julia> p
 PA1|1 PB3|4 ZE1 ZE2

julia> p.word
3-element Array{Tuple{Int64,Array{Operator,1}},1}:
 (1, [P1|1])
 (2, [P3|4])
 (5, [Z1, Z2])
```
It is assumed that:

1. parties are numbered starting from 1,
2. the party numbers are in strictly increasing order `p1 < p2 < p3 ...`, and
3. only parties that have at least one operator associated with them appear
   in the list.

This (mainly the definition of `Monomial` and the function
`Base.:*(x::Monomial, y::Monomial)`) is the part of the code that would have
to be modified if we wanted to support, e.g., operators acting on the
subystems of more than one party.
