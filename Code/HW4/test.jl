using LinearAlgebra

function spquad(dim, ord, bpt=nothing)
    if bpt !== nothing && size(bpt, 2) != dim
        throw(ArgumentError("dimension mismatch"))
    elseif bpt === nothing
        bpt = repeat([-1 1]', dims=1)  # Defaults to [-1, 1]^dim
    end

    len = diff(bpt, dims=1)
    mpt = mean(bpt, dims=1)

    if ord == 0
        x = mpt'
        w = prod(len)
    elseif ord != 0 && dim == 1
        x = (mpt .+ len .* cos.(π * 2 .^ (-ord) .* (0:2^ord)' ./ 2))'
        w = clencurt(2^ord + 1)' * len / 2
    else
        x0 = mpt
        w0 = prod(len)

        for j in 1:ord
            x, w = sparsegridnd(dim, j, mpt, len)
            I, = intersect(x, x0, dims=1)
            w[I] .+= w0
            x0, w0 = x, w
        end
    end

    return x, w
end

function sparsegridnd(n, ord, mpt, len)
    p(i) = unique((i .!= 0) .* cos.(π * 2 .^ (-i) .* (0:2^i)'))

    v = genindex(n, ord)
    vmax = 1 + maximum(v[:, n])

    P = uaf(j -> p(j), 0:vmax-1)
    Q = uaf(j -> diffweight(j), 0:vmax-1)

    m = size(v, 1)
    xw = getpts(P, Q, v, mpt, len)
    x = xw[:, 1:n]'
    w = prod(uaf(k -> xw[:, n+k], 1:length(P)), dims=1)

    x = round.(x, digits=15)

    _, ii, jj = unique(x, dims=2, return_index=true)
    x, w = x[:, ii], sparse(jj, 1:size(x, 2), 1, size(x, 2), size(x, 2)) * w

    return x, w
end

function diffweight(ii)
    if ii == 0
        return 2
    elseif ii == 1
        return [1; -2; 1] / 3
    else
        q(i) = clencurt(2^i + 1)'
        dw = q(ii)
        dw[1:2:end] .-= q(ii-1)
        return dw
    end
end

function clencurt(N1)
    if N1 == 1
        return [2]
    else
        N = N1 - 1
        c = zeros(N1)
        c[1:2:N1] .= 2 ./ [1; 1 .- (2:2:N).^2]
        f = real(ifft([c; reverse(c[2:N])]))
        return 2 * [f[1]; 2 * f[2:N]; f[N1]] / 2
    end
end

function genindex(n, L1, head=nothing)
    if n == 1
        return L1
    else
        return uaf(j -> genindex(n-1, L1-j, j), 0:L1)
    end
end

function getpts(P, Q, v, mpt, len)
    x = [mpt[k] .+ xk[:].*len[k] ./ 2 for (k, xk) in enumerate.(P[1 .+ v])]
    w = prod([wk[:] for wk in Q[1 .+ v]], dims=1)
    return hcat(x..., w)
end

function uaf(arg, dex, varargin=nothing)
    result = map(arg, dex)
    if nargin() > 2
        result = hcat(result...)
    end
    return result
end
