# Convert party numbers to strings and vice versa.

alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z']



"Return string representation of a party, e.g., `party_str(3) = \"C\"`."
function party_str(p::Integer)
    base = length(alphabet)
    chars = Array{Char,1}()

    while (p > 0)
        p -= 1
        push!(chars, alphabet[1 + p % base])
        p = div(p, base)
    end

    return String(reverse!(chars))
end

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
party_num(n::Integer) = n

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
