# bff-npa

Code to do NPA needed to use the Brown-Fawzi-Fawzi method.

At the moment: `bff.jl` contains code in development, `ops.jl` is the older
code that only does NPA with projectors and shouldn't be used (except to
copy/adapt some of the code in it and move to `bff.jl`).

Use `bff.jl` like this:
```julia
julia> using Printf

julia> using Base.Iterators

julia> include("bff.jl")
```


## Basic features

We can do arithmetic with and take conjugates of some different types of
operators that we associate to different parties. At the moment:
- dichotomic,
- fourier,
- projector,
- unitary,
- zbff (Brown-Fawzi-Fawzi operators).
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

julia> ZE = zbff(5, 1:2)
2-element Array{Monomial,1}:
 ZE1
 ZE2

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
projector(party, output, input)
unitary(party, index, conj=false)
zbff(party, index, conj=false)
```
In all of these:
- Party numbers start from 1.
- The parameters called `input`, `output`, and `index` can be either integers
  or arrays or ranges of integers.
- The parameter `conj` is optional and defaults to `false` if it is omitted.

Some examples:
```julia
julia> projector(1, 1:2, 1:2)
2×2 Array{Monomial,2}:
 PA1|1  PA1|2
 PA2|1  PA2|2

julia> julia> zbff(1, 1:3)
3-element Array{Monomial,1}:
 ZA1
 ZA2
 ZA3
```

Note that there are no special relations (at least, at the moment) between
these different types of operators, so you shouldn't mix, e.g., projectors
and dichotomic operators unless you consider them to be unrelated to each
other:
```julia
julia> dichotomic(1, 1) * projector(1, 1, 1)
A1 PA1|1

julia> dichotomic(1, 1) - (2*projector(1, 1, 1) - Id)
Id + A1 - 2 PA1|1
```

Note that monomials and polynomials are different types (although they are
occasionally printed the same), and it is possible to loop over the monomials
and (nonzero) coefficients in a polynomial:
```julia
julia> projector(1,1,1)
PA1|1

julia> typeof(projector(1,1,1))
Monomial

julia> 1*projector(1,1,1)
PA1|1

julia> typeof(1*projector(1,1,1))
Polynomial

julia> S
A1 B1 + A1 B2 + A2 B1 - A2 B2

julia> typeof(S)
Polynomial

julia> for (m, c) in S
           @printf "%s  =>  %2d\n" m c
       end
A2 B1  =>   1
A1 B2  =>   1
A2 B2  =>  -1
A1 B1  =>   1

julia> for (m, c) in sort(S)
           @printf "%s  =>  %2d\n" m c
       end
A1 B1  =>   1
A1 B2  =>   1
A2 B1  =>   1
A2 B2  =>  -1
```


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

for (i, x) in enumerate(ops2)
    for j in i:length(ops2)
        y = ops2[j]
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
inequalities) can also be used to compare polynomials.


## Internal details

The way a list of operators are joined to multiply them is determined at the
moment by a function `join_ops` near the beginning of the file `bff.jl`. It
is not super general at the moment but is general enough to handle the
different types of operators already defined.

Associating operators in groups to parties is handled by the `Monomial`
type. At the moment it just contains a list
```julia
word = [(p1, ops1), (p2, ops2), ...]
```
of parties `p1`, `p2`, etc. and lists of operators `ops1`, `ops2`,
etc. associated to those parties. For example,
```julia
julia> R
PA1|1 PB2|2 ZE1 ZE2

julia> R.word
3-element Array{Tuple{Integer,Array{Operator,1}},1}:
 (1, [P1|1])
 (2, [P2|2])
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
