# using Printf

alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z']



function party2string(p::Int)
    base = length(alphabet)
    chars = Array{Char,1}()

    while (p > 0)
        p -= 1
        push!(chars, alphabet[1 + p % base])
        p = div(p, base)
    end

    return String(reverse!(chars))
end



struct Monomial
    projectors::Array{Tuple{Int,Array{Tuple{Int,Int},1}},1}
end

Id = Monomial([])

function Base.show(io::IO, x::Monomial)
    if length(x.projectors) == 0
        print(io, " Id")
    else
        for (p, oips) in x.projectors
            for (output, input) in oips
                @printf io " %s%d|%d" party2string(p) output input
            end
        end
    end
end

Base.hash(x::Monomial) = hash(x.projectors)



struct Polynomial
    terms::Dict{Monomial,Number}
end

Polynomial() = Polynomial(Dict{Monomial,Number}())

Polynomial(x::Number) = Polynomial((x != 0) ? Dict(Id => x) : Dict())

Polynomial(x::Monomial) = Polynomial(Dict(x => 1))

Polynomial(x::Polynomial) = x

Polynomial(x::Base.Generator) = Polynomial(Dict(x))



function Base.copy(x::Polynomial)
    return Polynomial(copy(x.terms))
end

Base.iterate(x::Polynomial) = iterate(x.terms)
Base.iterate(x::Polynomial, state) = iterate(x.terms, state)

function Base.show(io::IO, p::Polynomial)
    terms = p.terms

    if isempty(terms)
        print(io, " 0")
    else
        for (m, c) in terms
            print(io, " + (", c, ")")
            show(io, m)
        end
    end
end


IndexRange = Union{UnitRange{Int},
                   StepRange{Int,Int},
                   Array{Int}}



function projector(party, output::Int, input::Int)
    return Monomial([(party, [(output, input)])])
end

function projector(party, output::IndexRange, input::Int)
    return [projector(party, o, input) for o in output]
end

function projector(party, output::IndexRange, input::IndexRange)
    return [projector(party, o, i) for o in output, i in input]
end


Base.:(==)(x::Number, y::Monomial) = (x == 1) && isempty(y.projectors)

Base.:(==)(x::Monomial, y::Number) = (y == x)

Base.:(==)(x::Monomial, y::Monomial) = (x.projectors == y.projectors)

function Base.:(==)(x::Number, y::Polynomial)
    return y.terms == Dict(Id => x)
end

Base.:(==)(x::Polynomial, y::Number) = (y == x)

function Base.:(==)(x::Number, y::Polynomial)
    return y.terms == Dict()
end

Base.:(==)(x::Polynomial, y::Polynomial) = (x.terms == y.terms)



function Base.getindex(x::Polynomial, y::Monomial)
    return (y in keys(x.terms)) ? x.terms[y] : 0
end

function Base.setindex!(x::Polynomial, y::Number, z::Monomial)
    if y == 0
        delete!(x.terms, z)
    else
        x.terms[z] = y
    end
end



function Base.:+(x::Number, y::Monomial)
    z = Polynomial(y)
    z += x
    return z
end

Base.:+(x::Monomial, y::Number) = y + x

function Base.:+(x::Monomial, y::Monomial)
    return Polynomial((x != y) ? Dict(x => 1, y => 1) : Dict(x => 2))
end

function Base.:+(x::Number, y::Polynomial)
    z = copy(y)
    z[Id] += x
    return z
end

Base.:+(x::Polynomial, y::Number) = y + x

function Base.:+(x::Monomial, y::Polynomial)
    z = copy(y)
    z[x] += 1
    return z
end

Base.:+(x::Polynomial, y::Monomial) = y + x

function Base.:+(x::Polynomial, y::Polynomial)
    z = copy(x)

    for (m, c) in y
        z[m] += c
    end

    return z
end



function Base.:-(x::Number, y::Monomial)
    z = Polynomial(Dict(y => -1))
    z += 1
    return z
end

function Base.:-(x::Monomial, y::Number)
    z = Polynomial(x)
    z -= y
    return z
end

function Base.:-(x::Monomial, y::Monomial)
    return Polynomial((x != y) ? Dict(x => 1, y => -1) : Dict())
end

function Base.:-(x::Number, y::Polynomial)
    z = Polynomial((m, -c) for (m, c) in y)
    z[Id] += 1
    return z
end

function Base.:-(x::Polynomial, y::Number)
    z = copy(x)
    z[Id] -= y
    return z
end

function Base.:-(x::Monomial, y::Polynomial)
    z = Polynomial((m, -c) for (m, c) in y)
    z[x] += 1
    return z
end

function Base.:-(x::Polynomial, y::Monomial)
    z = copy(x)
    z[y] -= 1
    return z
end

function Base.:-(x::Polynomial, y::Polynomial)
    z = copy(x)

    for (m, c) in y
        z[m] -= c
    end

    return z
end



function Base.:*(x::Number, y::Monomial)
    return (x != 0) ? Polynomial(Dict(y => x)) : 0
end

function Base.:*(x::Monomial, y::Number)
    return (y != 0) ? Polynomial(Dict(x => y)) : 0
end

function Base.:*(x::Monomial, y::Monomial)
    M = length(x.projectors)

    if M == 0
        return y
    end

    N = length(y.projectors)

    if N == 0
        return x
    end

    j = 1
    k = 1

    projectors = Array{Tuple{Int,Array{Tuple{Int,Int},1}},1}()

    while j <= N && k <= M
        (px, oipsx) = x.projectors[j]
        (py, oipsy) = y.projectors[k]

        if px < py
            push!(projectors, x.projectors[j])
            j += 1
        elseif py < px
            push!(projectors, y.projectors[k])
            k += 1
        else
            (ox, ix) = oipsx[end]
            (oy, iy) = oipsy[1]

            if ix == iy
                if ox == oy
                    push!(projectors, (px, vcat(oipsx, oipsy[2:end])))
                else
                    return 0
                end
            else
                push!(projectors, (px, vcat(oipsx, oipsy)))
            end

            j += 1
            k += 1
        end
    end

    append!(projectors, x.projectors[j:end])
    append!(projectors, y.projectors[k:end])

    return Monomial(projectors)
end

function Base.:*(x::Number, y::Polynomial)
    return (x != 0) ? Polynomial((m, x*c) for (m, c) in y) : 0
end

Base.:*(x::Polynomial, y::Number) = y * x

function Base.:*(x::Monomial, y::Polynomial)
    z = Polynomial()

    for (my, c) in y
        if (m = x*my) != 0
            z[m] += c
        end
    end

    return z
end

function Base.:*(x::Polynomial, y::Monomial)
    z = Polynomial()

    for (mx, c) in x
        if (m = mx*y) != 0
            z[m] += c
        end
    end

    return z
end

function Base.:*(x::Polynomial, y::Polynomial)
    z = Polynomial()

    for (mx, cx) in x
        for (my, cy) in y
            if (m = mx * my) != 0
                z += (cx*cy) * m
            end
        end
    end

    return z
end



Base.:/(x::Monomial, y::Number) = Polynomial(Dict(x => 1/y))

function Base.:/(x::Polynomial, y::Number)
    return Polynomial((m, c/y) for (m, c) in x)
end


function Base.conj(x::Monomial)
    return Monomial([(p, reverse(oips)) for (p, oips) in x.projectors])
end

function Base.adjoint(x::Monomial)
    return Monomial([(p, reverse(oips)) for (p, oips) in x.projectors])
end

function Base.conj(x::Polynomial)
    return Polynomial((conj(m), conj(c)) for (m, c) in p)
end

function Base.adjoint(x::Polynomial)
    return Polynomial((adjoint(m), conj(c)) for (m, c) in p)
end
