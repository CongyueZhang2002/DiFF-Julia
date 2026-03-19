module MultiGKTable

using QuadGK

export integration_settings, 
       integrate, 
       get_value

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

# GK components

mutable struct integration_settings{FF,FP,FN,PR,FB,TP}

    f_factor::FF
    f_Pert::FP
    f_NP::FN

    product_rule::PR
    pert_node_type::Type{TP}
    f_bounds::FB   

    orders::Vector{Int}
    N_maxs::Vector{Int}
    rtols::Vector{Float64}
    atols::Vector{Float64}
end

struct interval{BK}

    layer::Int
    a::Float64
    b::Float64

    vars::Vector{Float64}
    nodes::Vector{Float64}
    blocks::Vector{BK}

    I::Float64
    err::Float64
end

@generated function get_block_type(::Val{layer}, ::Type{pert_node_type}) where {layer, pert_node_type}
    type = :pert_node_type
    for _ in 2:layer
        type = :(Vector{interval{$type}})
    end
    return type
end

@inline function get_worst_interval(intervals::Vector{<:interval})
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

@inline function create_interval(settings::integration_settings, 
    a:: Float64, b::Float64, vars::Vector{Float64}, weight::Float64, layer::Int)

    # GK grid components
    order = settings.orders[layer]
    nodes_norm, weights_norm, weights_G_norm = get_GK_grid_norm(order)
    m = (a + b) * 0.5
    h = (b - a) * 0.5
    N_nodes = length(nodes_norm)

    # initialize
    IK = 0.0
    IG = 0.0
    i_G = 1
    nodes = (layer == 1) ? Vector{Float64}(undef, N_nodes) : Float64[]
    block_type = get_block_type(Val(layer), settings.pert_node_type)
    blocks = Vector{block_type}(undef, N_nodes)

    # reuse vars_in vector
    N_vars = length(vars)
    vars_in = Vector{Float64}(undef, N_vars + 1)
    @inbounds copyto!(vars_in, 2, vars, 1, N_vars) # creates [?, vars[1], ...]

    @inbounds for i_node in 1:N_nodes

        #build GK grid
        x = m + h * nodes_norm[i_node]
        vars_in[1] = x
        weights_K = h * weights_norm[i_node]
        weight_in = weight*weights_K

        if layer == 1
            f_Pert_value = settings.f_Pert(vars_in)
            f_NP_value = settings.f_NP(vars_in)
            value = settings.product_rule(f_Pert_value, f_NP_value)
            block = weight_in .* f_Pert_value
            nodes[i_node] = x
        elseif layer == 2
            factor = settings.f_factor(vars_in)
            a_in = settings.f_bounds[layer-1][1](vars_in)
            b_in = settings.f_bounds[layer-1][2](vars_in)
            value_raw, block = adaptive_GK_integration(settings=settings, a=a_in, b=b_in, vars=copy(vars_in), weight=factor*weight_in, layer=layer-1)
            value = factor * value_raw
        else
            a_in = settings.f_bounds[layer-1][1](vars_in)
            b_in = settings.f_bounds[layer-1][2](vars_in)
            value, block = adaptive_GK_integration(settings=settings, a=a_in, b=b_in, vars=copy(vars_in), weight=weight_in, layer=layer-1)
        end

        IK += weights_K * value
        blocks[i_node] = block

        if iseven(i_node)
            weights_G = h * weights_G_norm[i_G]
            IG += weights_G * value
            i_G += 1
        end
    end

    err = abs(IK - IG)
    return interval{block_type}(layer, a, b, vars, nodes, blocks, IK, err)
end

@inline function split(settings::integration_settings, 
    iv::interval, vars::Vector{Float64}, weight::Float64, layer::Int)
    m = (iv.a + iv.b) * 0.5
    iv_1 = create_interval(settings, iv.a, m, vars, weight, layer)
    iv_2 = create_interval(settings, m, iv.b, vars, weight, layer)
    return iv_1, iv_2
end

function adaptive_GK_integration(; settings::integration_settings, 
    a::Float64, b::Float64, vars::Vector{Float64}, weight::Float64, layer::Int)

    iv_0 = create_interval(settings, a, b, vars, weight, layer)
    intervals = Vector{typeof(iv_0)}(undef, 1)
    intervals[1] = iv_0

    I_total    = iv_0.I
    aerr_total = iv_0.err
    rerr_total = aerr_total / max(abs(I_total), 1e-12)

    rtol  = settings.rtols[layer]
    atol  = settings.atols[layer]
    N_max = settings.N_maxs[layer]

    while (rerr_total > rtol) && (aerr_total > atol) && (length(intervals) < N_max)

        i = get_worst_interval(intervals)
        interval_worst = intervals[i]

        subinterval_1, subinterval_2 = split(settings, interval_worst, vars, weight, layer)
        intervals[i] = subinterval_1
        insert!(intervals, i+1, subinterval_2)

        I_total += (subinterval_1.I + subinterval_2.I - interval_worst.I)
        aerr_total += (subinterval_1.err + subinterval_2.err - interval_worst.err)
        rerr_total = aerr_total / max(abs(I_total), 1e-12)
    end

    return I_total, intervals
end

function integrate(settings::integration_settings{FF,FP,FN,PR,FB,TP}) where {FF,FP,FN,PR,FB,TP}

    layer = length(settings.f_bounds)
    a = settings.f_bounds[layer][1]
    b = settings.f_bounds[layer][2]

    I_total, intervals = adaptive_GK_integration(
        settings=settings, a=a, b=b, vars=Float64[], weight=1.0, layer=layer)

    return I_total, get_table(settings, intervals)
end

# return pert table

function get_table(settings::integration_settings{FF,FP,FN,PR,FB,TP},
                   intervals::Vector{<:interval}) where {FF,FP,FN,PR,FB,TP}

    intervals_flattened = [intervals]

    while intervals_flattened[1][1].layer > 1
        intervals_flattened = reduce(vcat, (iv.blocks for leaf in intervals_flattened for iv in leaf))
    end

    TableEntry = NamedTuple{
        (:vars, :nodes, :pert_vec),
        Tuple{Vector{Float64}, Vector{Float64}, Vector{TP}}
    }
    table = Vector{TableEntry}()
    sizehint!(table, length(intervals_flattened))

    for leaf in intervals_flattened
        vars = leaf[1].vars

        nodes    = Float64[]
        pert_vec = Vector{TP}()   

        for iv in leaf
            append!(nodes, iv.nodes)
            append!(pert_vec, iv.blocks)
        end

        push!(table, (vars = vars, nodes = nodes, pert_vec = pert_vec))
    end

    return table
end

function get_value(settings, table)
    value = 0.0
    for entry in table

        N_vars = length(entry.vars)
        vars_in = Vector{Float64}(undef, N_vars + 1)
        @inbounds copyto!(vars_in, 2, entry.vars, 1, N_vars)

        @inbounds for i in eachindex(entry.nodes)
            vars_in[1] = entry.nodes[i]          
            NP = settings.f_NP(vars_in)
            value += settings.product_rule(entry.pert_vec[i], NP)
        end
    end
    return value
end

# QuadGK for comparison

function integrate_quadgk(settings::integration_settings)

    layers = length(settings.f_bounds)

    bound_value(f, vars) = f isa Function ? f(vars) : f

    function integrate_layer(layer::Int, vars::Vector{Float64})
        f_lo, f_up = settings.f_bounds[layer]
        a = bound_value(f_lo, vars)
        b = bound_value(f_up, vars)

        vars_in(x) = vcat(x, vars)

        integrand = (
            if layer == 1
                x -> begin
                    v = vars_in(x)
                    settings.product_rule(settings.f_Pert(v), settings.f_NP(v))
                end
            elseif layer == 2
                x -> begin
                    v = vars_in(x)
                    settings.f_factor(v) * integrate_layer(layer - 1, v)
                end
            else
                x -> integrate_layer(layer - 1, vars_in(x))
            end
        )

        val, _ = quadgk(integrand, a, b; rtol=settings.rtols[layer], atol=settings.atols[layer])

        return val
    end

    return integrate_layer(layers, Float64[])
end

end # module MultiGKTable

# test

test = 0
if test == 0
    nothing
elseif test == 1

    using QuadGK
    using SpecialFunctions
    using .MultiGKTable

    Q  = 50.0
    qT = min(0.2 * Q, 20.0)

    f_Pert(b) = b[1] * besselj0(qT*b[1]) 
    f_NP(b) = (Q/b[1])^(-0.02*b[1]^2)

    product_rule(a, b) = a * b

    bmin = 0.001
    bmax = 5 + 20 * exp(-Q / 50)

    settings = MultiGKTable.integration_settings(
        1.0,
        f_Pert,
        f_NP,
        product_rule,
        Float64,                        # pert_node_type
        [(bmin, bmax)],     
        [7],                            # orders
        [8],                            # N_maxs
        [0.001],                        # rtols
        [1e-8],                         # atols
    )

    val_table, table = MultiGKTable.integrate(settings)
    val_reconstructed = MultiGKTable.get_value(settings, table)
    val_ref, _ = quadgk(b -> (b * besselj0(qT * b)) * ((Q / b)^(-0.02 * b^2)),
                        bmin, bmax; rtol=1e-12)

    println("adaptive_GK value           = ", val_table)
    println("table reconstruction         = ", val_reconstructed)
    println("QuadGK reference             = ", val_ref)

    denom = max(abs(val_ref), 1e-300)
    println("diff(table - adaptive_GK)    = ", abs(val_reconstructed - val_table) / denom)
    println("rel diff(adaptive_GK - ref)  = ", abs(val_table - val_ref) / denom)

    num_nodes = sum(length(entry.nodes) for entry in table)
    println("num nodes = ", num_nodes)

elseif test == 2

    f_Pert(v) = (y=v[1]; x=v[2]; sqrt(max(1 - x^2 - y^2, 0.0)))
    f_NP(v) = 2.0
    product_rule(a,b) = a*b

    f_lo_y(v) = -sqrt(max(1 - v[1]^2, 0.0))
    f_up_y(v) = +sqrt(max(1 - v[1]^2, 0.0))

    f_factor(v) = 1.0

    settings = MultiGKTable.integration_settings(
        f_factor, f_Pert, f_NP, 
        product_rule, Float64,
        [(f_lo_y, f_up_y),(-1.0, 1.0)],          
        [7, 7], 
        [3, 3],            
        [1e-25, 1e-5], 
        [1e-12, 1e-12] 
    )

    V_num, table = MultiGKTable.integrate(settings)
    V_true = 4pi/3
    println("relerr = ", abs(V_num - V_true)/abs(V_true))

    table_value = MultiGKTable.get_value(settings, table)
    println("taberr = ", abs(table_value - V_num)/abs(V_num))

    value_gk = MultiGKTable.integrate_quadgk(settings)
    println("gkerr = ", abs(value_gk - V_true)/abs(V_true))
end