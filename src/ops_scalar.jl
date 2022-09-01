#abstract type Scalar

struct Trace
    word::PartiedOpList
end

function Base.show(io::IO, x::Trace)
    print(io, "Tr[")
    show(io, x.word)
    print(io, ']')
end

Base.hash(x::Trace, h::UInt) = hash(x.word, h)

function ctrace(m::Monomial)
    coeff = 1

    pcops = [(p, trace(ops)) for (p, ops) in m.word]

    for (_, (c, _)) in pcops
        coeff = rmul(coeff, c)
    end

    word = [(p, ops) for (p, (_, ops)) in pcops]

    return (coeff, word)
end



struct Ket
    label
end

struct Bra
    label
end

Base.show(io::IO, x::Ket) = print(io, '|', x.label, '>')
Base.show(io::IO, x::Bra) = print(io, '<', x.label, '|')

Base.conj(x::Ket) = Bra(x.label)
Base.conj(x::Bra) = Ket(x.label)

Base.hash(x::Ket, h::UInt) = hash(x.label, h)
Base.hash(x::Bra, h::UInt) = hash(x.label, h)



struct State
    word::PartiedOpList
    ket::Ket
end

struct CState
    bra::Bra
    word::PartiedOpList
end

ket(label) = State(PartiedOpList[], Ket(label)) 
bra(label) = CState(Bra(label), PartiedOpList[])

if !@isdefined(Psi)
    const Psi = ket(:Ψ)
end

function Base.show(io::IO, x::State)
    show_nonempty_pol(io, x.word)
    show(io, x.ket)
end

function Base.show(io::IO, x::CState)
    show(io, x.bra)
    show_nonempty_pol(io, x.word)
end

Base.conj(x::State) = CState(conj(x.ket), conj_pol(x.word))
Base.conj(x::CState) = State(conj_pol(x.word), conj(x.bra))

Base.hash(x::State, h::UInt) = hash((x.word, x.ket), h)
Base.hash(x::CState, h::UInt) = hash((x.bra, x.word), h)



struct BraKet
    bra::Bra
    word::PartiedOpList
    ket::Ket
end

function Base.show(io::IO, x::BraKet)
    word = x.word

    if isempty(word)
        print(io, '<', x.bra.label)
    else
        show(io, x.bra)
    end

    show_nonempty_pol(io, word)
    show(io, x.ket)
end

Base.conj(x::BraKet) = BraKet(conj(x.ket), conj_pol(x.word), conj(x.bra))
Base.hash(x::BraKet, h::UInt) = hash((x.bra, x.word, x.ket), h)

function Base.:*(x::CState, y::State)
    (c, word) = join_pol(x.word, y.word)

    @assert c == 1

    return BraKet(x.bra, word, y.ket)
end

function Base.:*(m::Monomial, x::State)
    (c, word) = join_pol(m.word, x.word)

    @assert c == 1

    return State(word, x.ket)
end

function Base.:*(x::CState, m::Monomial)
    (c, word) = join_pol(x.word, m.word)

    @assert c == 1

    return CState(x.bra, word)
end
