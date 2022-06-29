function gaussradau(m)
    x, v = FastGaussQuadrature.gaussradau(m);
    t = 0.5*(1 .- x);
    w = 0.5*v;

    return (t, w)
end

const ln2 = log(2)

zbff_op(M, Z, Zc, t) = M*(Z + Zc + (1 - t)*Zc*Z) + t*Z*Zc

function HAE_simple(Ms, constraints, level, m)
    ts, ws = gaussradau(m)

    # Drop the first one (make this optional later).
    ts = ts[2:end]
    ws = ws[2:end]

    Zs = zbff(5, 1:length(Ms))
    Zcs = conj(Z)

    result = 0

    for (t, w) in zip(ts, ws)
        c = w / (t*ln2)
        O = sum(zbff_op(M, Z, Zc, t)
                for (M, Z, Zc) in zip(Ms, Zs, Zcs))
        obj = npa_mind(O, constraints, level)
        result += c*(1 + obj)               
    end

    return result
end
