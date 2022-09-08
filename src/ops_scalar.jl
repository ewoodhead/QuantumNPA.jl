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

Base.conj(x::Trace) = Trace(conj_pol(x.word))
Base.adjoint(x::Trace) = conj(x)

function ctrace(m::Monomial)
    coeff = 1

    pcops = [(p, trace(ops)) for (p, ops) in m.word]

    for (_, (c, _)) in pcops
        coeff = rmul(coeff, c)
    end

    word = [(p, ops) for (p, (_, ops)) in pcops]

    return (coeff, word)
end



abstract type Ket end
abstract type Bra end

Base.show(io::IO, x::Ket) = print(io, '|', x.label, '>')
Base.show(io::IO, x::Bra) = print(io, '<', x.label, '|')

Base.hash(x::Ket, h::UInt) = hash(x.label, h)
Base.hash(x::Bra, h::UInt) = hash(x.label, h)

Base.adjoint(x::Ket) = conj(x)

struct AKet <: Ket
    label
end

struct ABra <: Bra
    label
end

Base.conj(x::AKet) = ABra(x.label)
Base.conj(x::ABra) = AKet(x.label)

struct OKet <: Ket
    label
    gram::Dict
end

struct OBra <: Bra
    label
    gram::Dict
end

Base.conj(x::OKet) = OBra(x.label, x.gram)
Base.conj(x::OBra) = OKet(x.label, x.gram)

function orthonormal_states(labels)
    gram = Dict(l => 1 for l in labels)
    return [OKet(l, gram) for l in labels]
end

struct GKet <: Ket
    label
    gram::Dict
end

struct GBra <: Bra
    label
    gram::Dict
end

Base.conj(x::GKet) = GBra(x.label, x.gram)
Base.conj(x::GBra) = GKet(x.label, x.gram)



struct State
    word::PartiedOpList
    ket::Ket
end

struct CState
    bra::Bra
    word::PartiedOpList
end

state(label) = State(PartiedOpList[], AKet(label))
cstate(label) = CState(ABra(label), PartiedOpList[])

if !@isdefined(Psi)
    const Psi = AKet(:Ψ)
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

Base.adjoint(x::State) = conj(x)
Base.adjoint(x::CState) = conj(x)

Base.hash(x::State, h::UInt) = hash((x.word, x.ket), h)
Base.hash(x::CState, h::UInt) = hash((x.bra, x.word), h)



struct BraKet
    bra::Bra
    word::PartiedOpList
    ket::Ket
end

braket(x::Bra, y::Ket) = BraKet(x, PartiedOpList[], y)

Base.:*(x::Bra, y::Ket) = braket(x, y)

function Base.:*(x::OBra, y::OKet)
    gram = x.gram

    if (y.gram === gram)
        l = x.label

        if (y.label != l)
            return 0
        elseif haskey(gram, l)
            return gram[l]
        end
    end

    return braket(x, y)
end

function Base.:*(x::GBra, y::GKet)
    gram = x.gram

    if (y.gram === gram)
        pair = (x.label, y.label)

        if haskey(gram, pair)
            return gram[pair]
        end
    end

    return braket(x, y)
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
Base.adjoint(x::BraKet) = conj(x)
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
