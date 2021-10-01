# using FastGaussQuadrature

function gaussradauln(m)
    x, v = gaussradau(m);
    t = 0.5*(1 .- x);
    w = 0.5*v;

    function f(x)
        y = x - 1
        return dot(w, y ./ (t*y .+ 1))
    end

    return f
end
