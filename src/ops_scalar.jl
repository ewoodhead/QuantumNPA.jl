#abstract type Scalar

struct Trace
    word::PartiedOpList
end

function Base.show(io::IO, x::Trace)
    print(io, "Tr[")
    show_pol(io, x.word)
    print(io, ']')
end

function ctrace(m::Monomial)
    coeff = 1

    pcops = [(p, trace(ops)) for (p, ops) in m.word]

    for (_, (c, _)) in pcops
        coeff = rmul(coeff, c)
    end

    x = Trace([(p, ops) for (p, (_, ops)) in pcops])

    return (coeff == 1) ? x : (coeff, x)
end
