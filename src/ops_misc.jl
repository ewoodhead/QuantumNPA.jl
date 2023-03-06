PartyVec = Vector{Int}



# Convert party numbers to strings and vice versa.

alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z']

party_sep = '_'



"Return string representation of a party, e.g., `party_str(3) = \"C\"`."
function party_str(p::Integer)
    base = length(alphabet)
    chars = Vector{Char}()

    while (p > 0)
        p -= 1
        push!(chars, alphabet[1 + p % base])
        p = div(p, base)
    end

    return String(reverse!(chars))
end

party_str(p::PartyVec) = join(map(party_str, p), party_sep)

party_str(s::String) = s
party_str(s::Symbol) = string(s)


pos(c::Char) = findfirst(isequal(c), alphabet)

"Return integer representation of a party, e.g., `party_num(\"C\") = 3`."
function party_num(s::AbstractString)
    base = length(alphabet)
    party = 0

    for c in s
        party = pos(c) + base * party
    end

    return party
end

party_num(c::Char) = pos(c)
party_num(s::Symbol) = party_num(string(s))
party_num(n::Int) = n

"""
Check if arg is a valid party vector, i.e., vector of ints where:
1) The first number is strictly greater than zero.
2) Every number is strictly greater than the previous one.
"""
function party_vec_p(party::PartyVec)
    if isempty(party)
        return true
    end

    (p0, ps) = Iterators.peel(party)

    if p0 < 1
        return false
    end

    for p in ps
        if p <= p0
            return false
        end

        p = p0
    end

    return true
end

function party_vec(s::AbstractString)
    party = map(party_num, split(s, party_sep))
    @assert party_vec_p(party)
    return party
end

party_vec(c::Char) = [party_num(c)]
party_vec(s::Symbol) = party(string(s))

function party_vec(n::Int)
    @assert n > 0
    return Int[n]
end

function party_vec(party::PartyVec)
    @assert party_vec_p(party)
    return party
end


"Split string into party and rest, e.g. \"AB1\" -> (\"AB\", \"1\")."
function split_party(s::String)
    k = findlast(in(alphabet), s)
    return (s[1:k], s[k+1:end])
end

split_party(s::Symbol) = split_party(string(s))



# Helper functions to print polynomials.

num2str(x::Real) = "$x"

function num2str(x::Rational)
    a, b = numerator(x), denominator(x)

    return (b != 1) ? "$a/$b" : "$a"
end

num2str(n::Bool) = (n ? "1" : "-1")

function csgn(x::Real, p::String = "+", m::String = "-")
    return (x >= 0) ? p : m
end

function sgnnum(x::Number, p::String = "+", m::String = "-")
    xr = real(x)
    xi = imag(x)
    
    if xi == 0
        return (csgn(xr, p, m), num2str(abs(xr)))
    elseif xr == 0
        return (csgn(xi, p, m), "$(num2str(abs(xi)))im")
    else
        xis = num2str(abs(xi))

        if xr >= 0
            xrs = num2str(xr)
            s = csgn(xi)
            
            return (p, "($xrs $s $(xis)im)")
        else
            xrs = num2str(-xr)
            s = csgn(-xi)

            return (m, "($xrs $s $(xis)im)")
        end
    end
end

function firstcoeff2string(x::Number)
    if x == 1
        return ""
    elseif x == -1
        return "-"
    else
        (s, xs) = sgnnum(x, "", "-")
        return "$s$xs "
    end
end

function coeff2string(x::Number)
    if x == 1
        return " + "
    elseif x == -1
        return " - "
    else
        (s, xs) = sgnnum(x)
        return " $s $xs "
    end
end
