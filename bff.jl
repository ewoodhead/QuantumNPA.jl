#using Printf

abstract type Operator end

# Default multiplication rule for operators, which can be specialised.
#
# Assumption for the moment: multiplication returns one of the three
# following things:
#   1) A pair (x, y) of operators.
#   2) A single operator x.
#   3) A number.
Base.:*(x::Operator, y::Operator) = (x, y)
Base.show(io::IO, o::Operator) = print_op(io, o)


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



struct Projector <: Operator
    output::Int
    input::Int
end

function print_op(io::IO, p::Projector, party=nothing)
    if party === nothing
        @printf io "P%d|%d" p.output p.input
    else
        @printf io "P%s%d|%d" party2string(party) p.output p.input
    end
end

function Base.:*(p::Projector, q::Projector)
    if p.input == q.input
        if p.output == q.output
            return p
        else
            return 0
        end
    else
        return (p, q)
    end
end

Base.conj(p::Projector) = p

struct Monomial
    coeff::Number
    word::Array{Tuple{Int64,Array{Operator,1}},1}
end

Id = Monomial(1, [])

function Base.show(io::IO, m::Monomial)
    print(io, m.coeff)
    if isempty(m.word)
        print(io, " Id")
    else
        for (party, ops) in m.word
            for o in ops
                print(io, " ")
                print_op(io, o, party)
            end
        end
    end
end

function Base.conj(m::Monomial)
    return Monomial(conj(m.coeff),
                    [(party, reverse!([conj(op) for op in ops]))
                     for (party, ops) in p.word])
end

function Base.:*(x::Number, y::Monomial)
    coeff = x * y.coeff
    return (coeff != 0) ? Monomial(coeff, y.word) : 0
end

Base.:*(x::Monomial, y::Number) = y * x

function Base.:*(x::Monomial, y::Monomial)
    coeff = x.coeff * y.coeff

    
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


function projector(party, output, input)
    return Monomial(1, [(party, [Projector(output, input)])])
end

