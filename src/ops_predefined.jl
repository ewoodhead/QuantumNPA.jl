@operator Dichotomic(input::Integer) ("$party$input", "/$input")

function Base.:*(x::Dichotomic, y::Dichotomic)
    return (x.input == y.input) ? (1, []) : (1, [x, y])
end

function parse_dichotomic(expr)
    if expr isa Symbol
        (party, input) = split_party(expr)
        return :($(esc(expr)) = dichotomic($party, $(parse(Int, input))))
    else
        name = expr.args[1]
        party = string(name)
        range = expr.args[2]
        return :($(esc(name)) = dichotomic($party, $(esc(range))))
    end
end

macro dichotomic(vars...)
    exprs = map(parse_dichotomic, vars)
    return quote $(exprs...) end
end



@operator(Fourier(input::Integer, power::Integer, d::Integer),
          "$party$input^$power",
          false,
          (d, input, power))

function Base.:*(x::Fourier, y::Fourier)
    input = x.input
    d = x.d

    if (y.d != d) || (y.input != input)
        return (1, [x, y])
    else
        p = (x.power + y.power) % d
        return (1, ((p != 0) ? [Fourier(input, p, d)] : []))
    end
end

Base.conj(x::Fourier) = Fourier(x.input, x.d - x.power, x.d)

function fourier(party, input::Integer, power::Integer, d::Integer)
    @assert d > 0
    p = power % d
    return (p != 0) ? Monomial(party, Fourier(input, p, d)) : Id
end

function fourier(party, input::IndexRange, power::Integer, d::Integer)
    return [fourier(party, z, power, d) for z in input]
end

function fourier(party, input::Integer, power::IndexRange, d::Integer)
    return [fourier(party, input, power, d) for p in power]
end

function fourier(party, input::IndexRange, power::IndexRange, d::Integer)
    return [fourier(party, z, p, d) for z in input, p in power]
end



@operator(Projector(output::Integer, input::Integer),
          "P$party$output|$input",
          (true, false))

function Base.:*(p::Projector, q::Projector)
    if p.input == q.input
        if p.output == q.output
            return (1, [p])
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [p, q])
    end
end

function projector(party, output::IndexRange, input::Integer;
                   full::Bool=false)
    if full
        outputs = output[1:end-1]
        ps = [projector(party, o, input) for o in outputs]
        pl = Polynomial(Id)

        for p in ps
            pl[p] = -1
        end

        return vcat(ps, [pl])
    else
        return [projector(party, o, input) for o in output]
    end
end

function projector(party, output::Integer, input::IndexRange)
    return [projector(party, output, i) for i in input]
end

function projector(party, output::IndexRange, input::IndexRange;
                   full::Bool=false)
    if full
        outputs = output[1:end-1]
        ps = [projector(party, o, i) for o in outputs, i in input]
        pl = reshape([Polynomial(Id) for _ in input], 1, length(input))

        for i in input
            for p in ps[:,i]
                pl[i][p] = -1
            end
        end

        return vcat(ps, pl)
    else
        return [projector(party, o, i) for o in output, i in input]
    end
end



@operator(KetBra(outputl::Integer, outputr::Integer, input::Integer),
          "|$outputl><$outputr|$input$party",
          (false, true))

function Base.:*(p::Projector, k::KetBra)
    if p.input == k.input
        if p.output == k.outputl
            return (1, [k])
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [p, k])
    end
end

function Base.:*(k::KetBra, p::Projector)
    if k.input == p.input
        if k.outputr == p.output
            return (1, [k])
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [k, p])
    end
end

function Base.:*(k::KetBra, l::KetBra)
    input = k.input

    if input == l.input
        kor, lol = k.outputr, l.outputl

        if kor == lol
            kol, lor = k.outputl, l.outputr

            if kol != lor
                return (1, [KetBra(kol, lor, input)])
            else
                return (1, [Projector(kol, input)])
            end
        else
            return (0, Array{Projector,1}())
        end
    else
        return (1, [k, p])
    end
end

Base.conj(k::KetBra) = KetBra(k.outputr, k.outputl, k.input)

function ketbra(party, outputl::Integer, outputr::Integer, input::Integer)
    op = ((outputl != outputr) ?
          KetBra(outputl, outputr, input) :
          Projector(outputl, input))
    return Monomial(party, op)
end



@operator Unitary(index::Integer, conj::Bool) "U$party$conj$index"

function Base.:*(u::Unitary, v::Unitary)
    if (u.index != v.index) || (u.conj == v.conj)
        return (1, [u, v])
    else
        return (1, [])
    end
end



@operator Generic(index::Integer, conj::Bool) "Z$party$conj$index"



@operator Hermitian(index::Integer) "H$party$index"



@operator(Pauli(k::Int, label::Char),
          "$party$label",
          false)

sigmaX = Pauli(1, 'X')
sigmaY = Pauli(2, 'Y')
sigmaZ = Pauli(3, 'Z')
pauli_operators = [sigmaX, sigmaY, sigmaZ]

function pauli(party, label::Char)
    k = findfirst(isequal(label), ['X', 'Y', 'Z'])

    if isnothing(k)
        @error "Pauli label must be 'X', 'Y', or 'Z'."
    end

    return Monomial(party, pauli_operators[k])
end

pauli_mul_table = [ (   1,    []   )   ( 1im, [sigmaZ])   (-1im, [sigmaY]);
                    (-1im, [sigmaZ])   (   1,    []   )   ( 1im, [sigmaX]);
                    ( 1im, [sigmaY])   (-1im, [sigmaX])   (   1,    []   )  ]

Base.:*(x::Pauli, y::Pauli) = pauli_mul_table[x.k, y.k]

function parse_pauli(var::Symbol)
    varname = string(var)

    party = varname[1:end-1]
    label = varname[end]

    return :($(esc(var)) = pauli($party, $label))
end

macro pauli(vars...)
    exprs = map(parse_pauli, vars)
    return quote $(exprs...) end
end
