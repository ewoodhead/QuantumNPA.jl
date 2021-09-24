# bff-npa

Code to do NPA needed to use the Brown-Fawzi-Fawzi method.

At the moment: `bff.jl` contains code in development, `ops.jl` is the older
code that only does NPA with projectors and shouldn't be used (except to
copy/adapt some of the code in it and move to `bff.jl`).

Use `bff.jl` like this:
```
using Printf
include("bff.jl")
```

At the moment it can do multiplications with projectors and those Z operators
in the Brown-Fawzi-Fawzi paper, e.g.,
```
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

julia> bffz(1,3)
 ZA3

julia> bffz(1,3,true)
 ZA*3

julia> conj(bffz(1,3))
 ZA*3

julia> conj(bffz(1,3,true))
 ZA3

julia> p = projector(1,1,1)*projector(2,3,4)*bffz(5,1)*bffz(5,2)
 PA1|1 PB3|4 ZE1 ZE2

julia> conj(p)
 PA1|1 PB3|4 ZE*2 ZE*1

julia> p*p
 PA1|1 PB3|4 ZE1 ZE2 ZE1 ZE2

julia> conj(p)*p
 PA1|1 PB3|4 ZE*2 ZE*1 ZE1 ZE2

julia> p*conj(p)
 PA1|1 PB3|4 ZE1 ZE2 ZE*2 ZE*1
```

Syntax is
```
projector(party, output, input)
bffz(party, index)
bffz(party, index, conj)
```
Party numbers start from 1.

The way a list of operators are joined to multiply them is determined at the
moment by a function `join_ops` near the beginning of the file. This is what
would need to be generalised if we wanted to support multiplication of
slightly more general types of operators. For example, if we introduced a
unitary type with `U* U = 1`, the function `join_ops` as it is at the moment
would simplify `U* U* U U` to `U* U` but not all the way to `1`.

Associating operators in groups to parties is handled by the `Monomial`
type. At the moment it just contains a list
```word = [(p1, ops1), (p2, ops2), ...]```
of parties `p1`, `p2`, etc. and lists of operators `ops1`, `ops2`,
etc. associated to those parties. It is assumed that 1) parties are numbered
starting from 1, 2) the party numbers are in strictly increasing order `p1 <
p2 < p3 ...`, and 3) only parties that have at least one operator associated
with them appear in the list. This (mainly the definition of `Monomial` and
the function `Base.:*(x::Monomial, y::Monomial)`) is the part of the code
that would have to be modified if we wanted to support, e.g., operators
acting on the subystems of more than one party.
