using DataStructures # ] add DataStructures
using LinearAlgebra

# Gauss–Kronrod(15/7) constants on [-1, 1]

const _xk = [
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.0
]
const _wk = [
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714
]
const _xg = [
    0.949107912342758524526189684047851,
    0.741531185599394439863864773280788,
    0.405845151377397166906606412076961,
    0.0
]
const _wg = [
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327
]

# One GK15 panel on [a,b]
# Returns: (Ik, Ig, err, nodes::Vector{Float64}, weights::Vector{Float64}, neval::Int=15)

function _gk15_panel(f, a::Float64, b::Float64)
    c = (a + b)/2
    h = (b - a)/2
    nξ = length(_xk) - 1  # number of positive Kronrod abscissae

    # evaluate center first to infer the accumulation type
    f0 = f(c)
    T0 = typeof(f0)

    # storage for +/- evaluations in the integrand's return type
    fpos = Vector{T0}(undef, nξ)
    fneg = Vector{T0}(undef, nξ)
    @inbounds for i in 1:nξ
        ξ = _xk[i]
        fpos[i] = f(c + h*ξ)
        fneg[i] = f(c - h*ξ)
    end

    # helper to create a zero of the same "shape" as f0
    zero_like(x) = x isa Number ? zero(x) : (zero(eltype(x)) .* x)

    # Kronrod (15-pt) estimate in the integrand's type
    Ik = zero_like(f0)
    @inbounds for i in 1:nξ
        Ik += _wk[i] * (fpos[i] + fneg[i])
    end
    Ik = h * (Ik + _wk[end]*f0)

    # Gauss (7-pt) estimate: Kronrod contains Gauss nodes at indices 2,4,6 and 0
    Ig = h * ( _wg[1]*(fpos[2] + fneg[2]) +
               _wg[2]*(fpos[4] + fneg[4]) +
               _wg[3]*(fpos[6] + fneg[6]) +
               _wg[4]*f0 )

    # scalar error estimate (norm for vector-valued)
    err = (T0 <: Number) ? abs(Ik - Ig) : LinearAlgebra.norm(Ik - Ig)

    # Emit local nodes/weights mapped to [a,b] (still Float64)
    m = 2*nξ + 1
    nodes   = Vector{Float64}(undef, m)
    weights = Vector{Float64}(undef, m)
    k = 1
    @inbounds for i in 1:nξ
        ξ = _xk[i]; w = _wk[i]
        nodes[k]   = c + h*ξ;  weights[k] = h*w; k += 1
        nodes[k]   = c - h*ξ;  weights[k] = h*w; k += 1
    end
    nodes[k] = c; weights[k] = h*_wk[end]

    return Ik, Ig, err, nodes, weights, 15
end

# Adaptive driver (stack-based binary bisection)

function quadgk_nodes(f, a::Real, b::Real; rtol=1e-8, atol=0.0, maxeval::Int=10^7, sort_nodes::Bool=true, max_nodes::Union{Nothing,Int}=nothing, min_leaves::Int=1)  # (defaults to old behavior)

    a1, b1 = float(a), float(b)
    (isfinite(a1) && isfinite(b1) && a1 ≠ b1) || throw(ArgumentError("quadgk_nodes requires a finite, non-degenerate [a,b]."))

    Ik, Ig, err, nloc, wloc, ncall = _gk15_panel(f, a1, b1)
    neval = ncall

    pq = PriorityQueue{NTuple{7,Any}, Float64}()
    enqueue!(pq, (a1, b1, Ik, Ig, err, nloc, wloc), -err)

    # each leaf will contribute 15 nodes (GK15)
    leaf_limit = isnothing(max_nodes) ? typemax(Int) : max(1, fld(max_nodes, 15))
    safety = 0.5  # be stricter than Esum <= rtol·|I| for oscillatory cases

    while true
        # build scalar and component-wise tests
        Esum = 0.0

        firstkey = first(keys(pq))
        TIk = typeof(firstkey[3])
        Iest = zero(TIk)      # signed estimate in the integrand’s type

        Ecomp = nothing       # per-component error (immutable-safe accumulator)
        for (ai, bi, Ikp, Igp, errp, nlocp, wlocp) in keys(pq)
            Esum += errp
            Iest += Ikp
            if !(Ikp isa Number)
                Δ = Ikp - Igp
                Ecomp = (Ecomp === nothing) ? abs.(Δ) : (Ecomp .+ abs.(Δ))
            end
        end

        # scalar tolerance from signed estimate (good for cancellations)
        Iproxy = Iest isa Number ? abs(Iest) : LinearAlgebra.norm(Iest)
        TOL    = atol + rtol * max(Iproxy, eps())

        # component-wise tolerance (only if vector-valued)
        ok_comp = true
        if Ecomp !== nothing
            TOLvec = (atol .+ rtol .* abs.(Iest))
            ok_comp = all(Ecomp .<= safety .* TOLvec)        # NEW: safety factor
        end

        leaf_count = length(keys(pq))

        # (1) stop if tolerances met AND we have at least min_leaves
        ok_scalar = (Esum ≤ safety * TOL)                    # NEW: safety factor
        if ((ok_scalar && ok_comp) && (leaf_count ≥ min_leaves)) || (neval ≥ maxeval)
            break
        end
        # (2) stop if we reached the leaf cap
        if leaf_count ≥ leaf_limit
            break
        end

        # refine worst panel
        key, _ = dequeue_pair!(pq)
        ai, bi, Ik_old, Ig_old, err_old, nloc_old, wloc_old = key

        ci = (ai + bi)/2
        IkL, IgL, errL, nlocL, wlocL, ncallL = _gk15_panel(f, ai, ci)
        IkR, IgR, errR, nlocR, wlocR, ncallR = _gk15_panel(f, ci, bi)
        neval += ncallL + ncallR

        enqueue!(pq, (ai, ci, IkL, IgL, errL, nlocL, wlocL), -errL)
        enqueue!(pq, (ci, bi, IkR, IgR, errR, nlocR, wlocR), -errR)
    end

    # collect accepted leaves
    firstkey = first(keys(pq))
    T = typeof(firstkey[3])
    I_total = zero(T); E_total = 0.0
    nodes = Float64[]; weights = Float64[]
    for (ai, bi, Ikp, Igp, errp, nlocp, wlocp) in keys(pq)
        I_total += Ikp;  E_total += errp
        append!(nodes, nlocp);  append!(weights, wlocp)
    end
    if sort_nodes
        p = sortperm(nodes); nodes = nodes[p]; weights = weights[p]
    end
    return I_total, E_total, nodes, weights, neval
end

### For QuadGK comptibility

# Finite [a,b]
function _quadgk_compat_finite(f, a::Real, b::Real; rtol::Real=1e-8, atol::Real=0.0, maxevals::Integer=10^7, order::Integer=15, norm = LinearAlgebra.norm, sort_nodes::Bool=false, max_nodes::Union{Nothing,Int}=nothing, bin_mode::Bool=false)
    # If the caller is doing a smooth "bin" and wants ≤30 nodes, use the bin builder.
    if bin_mode && !isnothing(max_nodes) && max_nodes ≤ 30
        I, E, _x, _w, _ne, _rel = bin_table_leq30(f, a, b;
            target_relerr = max(rtol, 0.01),   # aim at <= 1% or rtol, whichever larger
            rtol_hint     = rtol,
            verify        = false,             # flip true during bring-up if you like
            ref_rtol      = min(rtol*1e-3, 1e-12),
            sort_nodes    = sort_nodes)
        return I, E
    end

    # Otherwise use the general adaptive driver with an optional node budget
    I, E, _x, _w, _ne = quadgk_nodes(f, a, b;
                                     rtol=rtol, atol=atol,
                                     maxeval=maxevals,
                                     sort_nodes=sort_nodes,
                                     max_nodes=max_nodes)
    return I, E
end

function _quadgk_compat_nodes(f, a::Real, b::Real; rtol::Real=1e-8, atol::Real=0.0, maxevals::Integer=10^7, order::Integer=15, norm = LinearAlgebra.norm, sort_nodes::Bool=false, max_nodes::Union{Nothing,Int}=nothing, bin_mode::Bool=false)
    # If the caller is doing a smooth "bin" and wants ≤30 nodes, use the bin builder.
    if bin_mode && !isnothing(max_nodes) && max_nodes ≤ 30
        I, E, _x, _w, _ne, _rel = bin_table_leq30(f, a, b;
            target_relerr = max(rtol, 0.01),   # aim at <= 1% or rtol, whichever larger
            rtol_hint     = rtol,
            verify        = false,             # flip true during bring-up if you like
            ref_rtol      = min(rtol*1e-3, 1e-12),
            sort_nodes    = sort_nodes)
        return I, E
    end

    # Otherwise use the general adaptive driver with an optional node budget
    I, E, x, w, _ne = quadgk_nodes(f, a, b;
                                     rtol=rtol, atol=atol,
                                     maxeval=maxevals,
                                     sort_nodes=sort_nodes,
                                     max_nodes=max_nodes)
    return x, w
end

# tests
test = false
if test
    using QuadGK
    using SpecialFunctions

    Q = 50.0
    qT = min(0.2*Q, 20.0)

    f(b) = b * besselj0(qT*b) * (Q/b)^(-0.02*b^2)

    bmin, bmax = 0.001, 5+20*exp(-Q/50)

    val_table, err = _quadgk_compat_finite(f, bmin, bmax; rtol=1e-2, max_nodes=100)

    val_ref, err_ref = quadgk(f, bmin, bmax; rtol=1e-12)

    println("rel diff            = ", abs(val_table - val_ref) / abs(val_ref))
end