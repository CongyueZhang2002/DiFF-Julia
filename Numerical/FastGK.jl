module FastGK

using QuadGK

export integrate

# Cache for QuadGK.kronrod 

const GK_grids_norm = Dict{
    Int,
    NamedTuple{
        (:nodes_norm,:weights_K_norm,:weights_G_norm),
        Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    }
}()

@inline function get_GK_grid_norm(order::Int)

    if haskey(GK_grids_norm, order)
        return GK_grids_norm[order]
    else
        nodes_norm, weights_K_norm, weights_G_norm = QuadGK.kronrod(order, -1.0, 1.0)
        grid = (; 
            nodes_norm = nodes_norm, 
            weights_K_norm = weights_K_norm, 
            weights_G_norm = weights_G_norm
        )
        GK_grids_norm[order] = grid
        return grid
    end
end

@inline function get_GK_grid(a::Float64, b::Float64,order::Int)

    grid_norm = get_GK_grid_norm(order)

    m = (a + b) * 0.5
    h = (b - a) * 0.5

    nodes = m .+ h .* grid_norm.nodes_norm
    weights_K = h .* grid_norm.weights_K_norm
    weights_G = h .* grid_norm.weights_G_norm

    return (nodes, weights_K, weights_G)
end    

struct interval
    a::Float64
    b::Float64
    I::Float64
    err::Float64
end

# GK components

@inline function create_interval(f, a::Float64, b::Float64, order::Int)

    nodes, weights_K, weights_G = get_GK_grid(a, b, order)
    N_nodes = length(nodes)

    IK = 0.0
    IG = 0.0
    i_G = 1

    @inbounds for i_node in 1:N_nodes

        x = nodes[i_node]
        value = f(x)
        IK += weights_K[i_node] * value

        if iseven(i_node)
            IG += weights_G[i_G] * value
            i_G += 1
        end
    end

    err = abs(IK - IG)
    return interval(a, b, IK, err)
end

@inline function split(f, iv::interval, order::Int)
    m = (iv.a + iv.b) * 0.5
    iv_1 = create_interval(f, iv.a, m, order)
    iv_2 = create_interval(f, m, iv.b, order)
    return iv_1, iv_2
end

@inline function get_worst_interval(intervals::Vector{interval})
    i_worst = 1
    err_max = intervals[1].err
    @inbounds for i in 2:length(intervals)
        err = intervals[i].err
        if err > err_max
            err_max = err
            i_worst = i
        end
    end
    return i_worst
end

function integrate(; 
    f,
    a::Float64, b::Float64,
    rtol::Float64, atol::Float64=1e-12,
    N::Int=10, order::Int=7
    )

    iv_0 = create_interval(f, a, b, order)
    intervals = [iv_0]  

    I_total    = iv_0.I
    aerr_total = iv_0.err
    rerr_total = aerr_total / max(abs(I_total), 1e-12)

    while (rerr_total > rtol) && (aerr_total > atol) && (length(intervals) < N)

        i = get_worst_interval(intervals)
        interval_worst = intervals[i]

        subinterval_1, subinterval_2 = split(f, interval_worst, order)
        intervals[i] = subinterval_1
        insert!(intervals, i+1, subinterval_2)

        I_total += (subinterval_1.I + subinterval_2.I - interval_worst.I)
        aerr_total += (subinterval_1.err + subinterval_2.err - interval_worst.err)
        rerr_total = aerr_total / max(abs(I_total), 1e-12)
    end

    return I_total
end

end # module FastGK

test = false
if test
    using .FastGK
    using QuadGK
    using SpecialFunctions

    Q = 50.0
    qT = min(0.2*Q, 20.0)

    f(b) = b * besselj0(qT*b)*(Q/b)^(-0.02*b^2)

    bmin, bmax = 0.001, 5+20*exp(-Q/50)

    val = FastGK.integrate(
        f=f,
        a=bmin, b=bmax,
        rtol=1e-2
    )

    val_ref, err_ref = quadgk(f, bmin, bmax; rtol=1e-12)

    println("adaptive_GK value      = ", val)
    println("QuadGK reference             = ", val_ref)
    println("rel diff            = ", abs(val - val_ref) / abs(val_ref))
end
