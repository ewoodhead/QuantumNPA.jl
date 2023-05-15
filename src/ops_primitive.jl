# Abstract type representing primitive operators.
# Specific ones are defined in the file 
abstract type Operator end

# Default multiplication rule for operators, which can be specialised.
#
# Assumption for the moment: multiplication returns a pair
#   (c, [ops...])
# consisting of a coefficient and a possibly empty list of operators.
#
# By default we return c = 1 and a list of the same two operators given as
# inputs, i.e., we just concatenate the operators.
Base.:*(x::Operator, y::Operator) = (1, [x, y])

# This controls how lists of operators are multiplied.
# It is not very general at the moment.
# Assumption: inputs opsx and opsy both contain at least one element.
function join_ops(opsx::Vector{Operator}, opsy::Vector{Operator})
    j = length(opsx)
    k = 1
    K = 1 + length(opsy)
    c = 1

    while true
        opx, opy = opsx[j], opsy[k]
        (c1, op) = opx * opy

        if c1 == 0
            return (0, Operator[])
        end

        c = rmul(c, c1)
        j -= 1
        k += 1

        if (op != []) || (j == 0) || (k == K)
            ops = Vector{Operator}(vcat(opsx[1:j], op, opsy[k:end]))
            return (c, ops)
        end
    end
end

# Default equality test. This should be specialised for specific types of
# operators (e.g., projectors), so if this default one is called it means the
# arguments are not the same type and are therefore not equal.
Base.:(==)(x::Operator, y::Operator) = false

# Default order. This again should be specialised for operators of the same
# type so here we just see if the type names are ordered lexicographically.
function Base.isless(x::Operator, y::Operator)
    return isless(nameof(typeof(x)), nameof(typeof(y)))
end

# This needs to be redefined for operators that aren't Hermitian.
Base.conj(o::Operator) = o

Base.adjoint(o::Operator) = conj(o)

Base.show(io::IO, o::Operator) = print(io, string(o))



"Simplify (if possible) Tr[ops...]."
function intrace_reduce(ops::Vector{Operator})
    c = 1

    while length(ops) > 2
        x = ops[end]
        y = ops[1]
        zs = ops[2:end-1]

        (c1, xy) = x*y

        if length(xy) >= 2
            break
        end

        c = rmul(c, c1)
        ops = Vector{Operator}(vcat(xy, zs))
    end

    return (c, ops)
end

function opcycles(ops::Vector{Operator})
    return (Vector{Operator}(vcat(ops[j:end], ops[1:j-1]))
            for j in 1:length(ops))
end

function trace(ops::Vector{Operator})
    (c, ops) = intrace_reduce(ops)
    return (c, minimum(opcycles(ops)))
end



# Definition of the @operator macro to define operators.

IndexRange = Union{UnitRange{<:Integer},Array{<:Integer}}



getfields(expr) = (expr.head == :tuple) ? expr.args : [expr]
getfieldnames(fields) = map(argname, fields)

argname(arg::Symbol) = arg
argname(arg::Expr) = arg.args[1]



instance_fields(instance, names) = [:($instance.$f) for f in names]

function fmt_remove(fmt, s::Symbol)
    return Expr(fmt.head, filter(!isequal(s), fmt.args)...)
end

function parse_fmt(fmt)
    if fmt.head === :string
        return (fmt, fmt_remove(fmt, :party))
    else
        fmts = fmt.args
        return (fmts[1], fmts[2])
    end
end

function stringfdef(name, fmt)
    fieldnames = [f for f in fmt.args if f isa Symbol]

    args = if (:party in fieldnames)
               (:(x::$name), :(party::PartyVec))
           else
               (:(x::$name),)
           end

    function fix_special(f)
        if f === :conj
            return :((x.$f ? "*" : ""))
        elseif f === :party
            return :(party_str(party))
        else
            return :(x.$f)
        end
    end
    
    bindings = (:($f = $(fix_special(f))) for f in fieldnames)

    return :(function Base.string($(args...))
                 $(bindings...)
                 return $fmt
             end)
end

function string_fdefs(name, fmt)
    (fmt_party, fmt_noparty) = parse_fmt(fmt)
    return (stringfdef(name, fmt_party), stringfdef(name, fmt_noparty))
end



function conjfalse(field)
    return ((argname(field) === :conj) ? Expr(:kw, field, false) : field)
end

function conj_def(name, fieldnames)
    if :conj in fieldnames
        cxfields = replace(instance_fields(:x, fieldnames),
                           :(x.conj) => :(!x.conj))
        conjf = :( Base.conj(x::$name) = $name($(cxfields...)) )
    else
        conjf = nothing
    end

    return conjf
end



parse_ctoropt(s::Bool) = (s, s)
parse_ctoropt(x::Expr) = x.args

function int_args(fields)
    return filter(fields) do f
        !isa(f, Symbol) && (f.args[2] === :Integer)
    end
end

function replace_seq(collection, replacements)
    collection = copy(collection)

    for r in replacements
        replace!(collection, r)
    end

    return collection
end

chtype(arg::Expr, type::Symbol) = Expr(arg.head, arg.args[1], type)
chtype(args::Array{Expr,1}, type::Symbol) = [chtype(a, type) for a in args]

function chtypes(fields, args::Array{Expr,1}, type::Symbol)
    return replace_seq(fields, (a => chtype(a, type) for a in args))
end

function ctor_ranges(lcname, fields, fieldnames)
    ifields = int_args(fields)
    to_sub = drop(powerset(ifields), 1)

    function mkrange(sub)
        fnames = getfieldnames(sub)
        nfields = chtypes(fields, sub, :IndexRange)
        lfcall = :($lcname(party, $(fieldnames...)))
        lassgms = [:($a = $a) for a in fnames]
        comp = Expr(:comprehension, Expr(:generator, lfcall, lassgms...))
        
        return :(function $lcname(party, $(nfields...))
                     return $comp
                 end)
    end

    return map(mkrange, to_sub)
end

function constructor_defs(make_constructors, name, fields, fieldnames)
    lcname = Symbol(lowercase(string(name)))
    (mk_ctor, mk_crange) = parse_ctoropt(make_constructors)

    cargs = Any[conjfalse(f) for f in fields]

    if mk_ctor
        mctor = :(function $(esc(lcname))(party, $(cargs...))
                      return OpProduct(party, $name($(fieldnames...)))
                  end)
    else
        mctor = nothing
    end

    if mk_crange
        mcrange = ctor_ranges(esc(lcname), cargs, fieldnames)
    else
        mcrange = ()
    end

    return (mctor, mcrange)
end



"""
Define a new type of operator with a given name (e.g., Projector), fields,
and format string. In addition to generating the struct definition this also
generates method definitions for the following generic functions:

  * Base.hash,
  * Base.:(==),
  * Base.isless,
  * Base.string,

as well as a constructor method with the name in lowercase (e.g., projector)
that creates a monomial containing a single operator associated to a given
party.

If one of the fields is named conj a Base.conj method is also generated.
"""
macro operator(ctor::Expr, fmt, make_constructors=true, order=nothing)
    name = ctor.args[1]
    fields = ctor.args[2:end]

    fieldnames = getfieldnames(fields)

    order = (isnothing(order) ? reverse(fieldnames) : getfields(order))
    xfields = :($(instance_fields(:x, order)...),)
    yfields = :($(instance_fields(:y, order)...),)

    (strf_party, strf_noparty) = string_fdefs(name, fmt)

    conjf = conj_def(name, fieldnames)
    (mctor, mcrange) = constructor_defs(make_constructors,
                                        name,
                                        fields,
                                        fieldnames)

    methods = [strf_party, strf_noparty, conjf, mctor, mcrange...]
    methods = filter(!isnothing, methods)

    return quote
        struct $name <: Operator
            $(fields...)
        end
        Base.hash(x::$name, h::UInt) = hash($xfields, h)
        Base.:(==)(x::$name, y::$name) = ($xfields == $yfields)
        Base.isless(x::$name, y::$name) = ($xfields < $yfields)
        $(methods...)
    end
end
