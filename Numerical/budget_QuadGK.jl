using QuadGK

function quadgk_nodes_weights_budget(f, a, b; rtol=1e-2, order=5, N=typemax(Int))
    @inline function eval_leaf(a1::Float64, b1::Float64)
        x, w, wg = kronrod(order, a1, b1)
        QK = 0.0; QG = 0.0; gi = 1
        @inbounds for j in 1:length(x)
            y = f(x[j])
            QK += w[j] * y
            if iseven(j); QG += wg[gi] * y; gi += 1; end
        end
        (a=a1, b=b1, x=x, w=w, QK=QK, E=abs(QK - QG))
    end

    a0 = Float64(a); b0 = Float64(b)
    leaf0 = eval_leaf(a0, b0)
    leaves = [leaf0]
    nodes_per_leaf = 2*order + 1
    I = leaf0.QK; E = leaf0.E

    while true
        tol = rtol * abs(I)
        if (E <= tol) || (length(leaves) >= N)
            sort!(leaves, by = l -> l.a)
            nL = length(leaves)
            xs = Vector{Float64}(undef, nL * nodes_per_leaf)
            ws = Vector{Float64}(undef, nL * nodes_per_leaf)
            pos = 1
            @inbounds for l in leaves
                copyto!(xs, pos, l.x, 1, nodes_per_leaf)
                copyto!(ws, pos, l.w, 1, nodes_per_leaf)
                pos += nodes_per_leaf
            end
            return xs, ws, I
        end

        # worst leaf
        i_max = 1; e_max = @inbounds leaves[1].E
        @inbounds for i in 2:length(leaves)
            ei = leaves[i].E
            if ei > e_max; e_max = ei; i_max = i; end
        end

        leaf = leaves[i_max]
        m = (leaf.a + leaf.b) * 0.5
        left  = eval_leaf(leaf.a, m)
        right = eval_leaf(m, leaf.b)

        I += (left.QK + right.QK - leaf.QK)
        E += (left.E  + right.E  - leaf.E)

        @inbounds leaves[i_max] = left
        push!(leaves, right)
    end
end

function quadgk_nodes_weights_budget_N1(f, a, b; rtol=1e-2, order=7, N=typemax(Int))
    a1 = Float64(a); b1 = Float64(b)
    x, w, _ = kronrod(order, a1, b1)
    return x, w
end

function quadgk_nodes_weights_budget_N2(f, a, b; rtol=1e-2, order=7, N=typemax(Int))
    a1 = Float64(a); b1 = Float64(b)

    if N <= 1
        x, w, _ = kronrod(order, a1, b1)
        return x, w
    end

    m = (a1 + b1) * 0.5
    xL, wL, _ = kronrod(order, a1, m)
    xR, wR, _ = kronrod(order, m, b1)

    nL = length(xL); nR = length(xR)
    xs = Vector{Float64}(undef, nL + nR)
    ws = Vector{Float64}(undef, nL + nR)

    @inbounds begin
        copyto!(xs, 1,     xL, 1, nL); copyto!(xs, nL+1, xR, 1, nR)
        copyto!(ws, 1,     wL, 1, nL); copyto!(ws, nL+1, wR, 1, nR)
    end

    return xs, ws
end

function bin_integrate(f, a, b; rtol_exact=1e-2, rtol_budget=5e-2, order=5, N=3)

    if bin_integrator == "QuadGK"
        return quadgk(f, a, b; rtol=rtol_exact)[1]
    elseif bin_integrator == "BudgetedQuadGK"
        xs, ws, val = quadgk_nodes_weights_budget(f, a, b; rtol=rtol_budget, order=order, N=N)
        return val
        #return weighted_sum(f, xs, ws)
        #xs, ws = [a,b], [0.5,0.5] # trapezoidal rule as a placeholder
    else
        error("Unknown bin_integrator: $bin_integrator")
    end
end

# --- quick demo ---
#    f(x) = (x-1)^2*exp(-x^2)
#    xs, ws, val = quadgk_nodes_weights_budget(f, 0.0, 10.0; rtol=1e-12, N=3)
#    I = sum(ws .* f.(xs))
#    # sanity: compare to quadgk (no budget)
#    I_ref, _ = quadgk(f, 0.0, 100.0; rtol=1e-12, atol=0.0)
#    println("|I - I_ref| = ", abs(I - I_ref)/abs(I_ref))
#    println("N nodes = ", length(xs))