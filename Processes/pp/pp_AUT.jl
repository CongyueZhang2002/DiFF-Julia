####################################################################################################
# Transversely Polarized 
####################################################################################################

function pp_EEC_T_integrand(; x2::Real, pT::Real, θ::Real, η::Real, sqrts::Real, rep::Int)

    # Hadronic                          
    S = sqrts^2
    T = -sqrts * pT * exp(-η)
    U = -sqrts * pT * exp(η)

    denom = x2 * S + T
    denom <= 0 && return 0.0

    # Partonic
    x1 = -x2 * U / denom
    (x1 <= 0 || x1 >= 1 || x2 <= 0 || x2 >= 1) && return 0.0

    s = Float64(x1 * x2 * S)
    t = Float64(x1 * T)
    u = Float64(x2 * U)

    # Scales
    μ = Float64(pT)
    Q = pT * exp(η)
    z = (1 - cos(θ)) / 2

    # Functions
    H_T = H_T_func(s=s, t=t, u=u)

    h1_vec = get_h1(x1, μ, rep)
    f1_vec = get_f1(x2, μ)
    H1a_vec = H1a_EEC_vec(zchi=z, Q=Q, mu=μ, rep=rep)

    # Sum over flavors
    val = 0.0

    for a in 1:11
        for b in 1:11
            for c in 1:11
                channel = H_T_channel(a=Flavors[a], b=Flavors[b], c=Flavors[c])
                channel === :zero && continue

                val += h1_vec[a] * f1_vec[b] * H1a_vec[c] * H_T[channel]
            end
        end
    end

    αs = αs_func(μ)
    factor = sin(θ)/2 * αs^2/(2*pi*S) * 1/(x1 * x2) * 1/denom

    return - factor * val
end

function pp_EEC_T(; pT::Real, θ::Real, η::Real, sqrts::Real, rep::Int)

    S = sqrts^2
    T = -sqrts * pT * exp(-η)
    U = -sqrts * pT * exp(η)

    denom_min = S + U
    denom_min <= 0 && return 0.0
    x2_min = min(max(-T / denom_min, 0.0), 1.0)

    integrand(x2) = pp_EEC_T_integrand(x2=x2, pT=pT, θ=θ, η=η, sqrts=sqrts, rep=rep)

    result = FastGK.integrate(f = integrand, a = x2_min, b = 1.0, rtol = 1e-2, N = 4)
    return result
end

function pp_EEC_T_bin(; pT::Real, θ_lo::Real, θ_hi::Real, η_lo::Real, η_hi::Real, sqrts::Real, rep::Int)

    S = sqrts^2

    x2_lo(v) = begin
        η = v[2]
        T = -sqrts * pT * exp(-η)
        U = -sqrts * pT * exp(η)
        denom = S + U
        denom <= 0 && return 1.0
        return min(max(-T / denom, 0.0), 1.0)
    end
    x2_hi(v) = 1.0

    θ_lo_func(v) = Float64(θ_lo)
    θ_up_func(v) = Float64(θ_hi)

    f_factor(v) = 1.0
    f_NP(v) = 1.0
    f_Pert(v) = pp_EEC_T_integrand(
        x2 = v[1],
        θ = v[2],
        η = v[3],
        pT = pT,
        sqrts = sqrts,
        rep = rep,
    )
    product_rule(pert, NP) = pert * NP

    settings = MultiGKTable.integration_settings(
        f_factor,
        f_Pert,
        f_NP,
        product_rule,
        Float64,
        [(x2_lo, x2_hi), (θ_lo_func, θ_up_func), (Float64(η_lo), Float64(η_hi))],
        [7, 5, 5],
        [5, 3, 3],
        [1e-2, 5e-2, 5e-2],
        [1e-12, 1e-12, 1e-12],
    )

    result, _ = MultiGKTable.integrate(settings)
    return result
end

####################################################################################################
# Unpolarized
####################################################################################################

function pp_EEC_U_integrand(; x2::Real, pT::Real, θ::Real, η::Real, sqrts::Real, rep::Int)

    # Hadronic
    S = sqrts^2
    T = -sqrts * pT * exp(-η)
    U = -sqrts * pT * exp( η)

    denom = x2 * S + T
    denom <= 0 && return 0.0

    # Partonic
    x1 = -x2 * U / denom
    (x1 <= 0 || x1 >= 1 || x2 <= 0 || x2 >= 1) && return 0.0

    s = Float64(x1 * x2 * S)
    t = Float64(x1 * T)
    u = Float64(x2 * U)

    # Scales
    μ = Float64(pT)
    Q = pT * exp(η)
    z = (1 - cos(θ)) / 2

    # Functions
    H_U = H_U_func(s=s, t=t, u=u)

    f1_x1_vec = get_f1(x1, μ)
    f1_x2_vec = get_f1(x2, μ)
    D1_vec = D1_EEC_vec(zchi=z, Q=Q, mu=μ, rep=rep)

    # Sum over flavors
    val = 0.0

    for a in 1:11
        for b in 1:11
            for c in 1:11
                channel = H_U_channel(a=Flavors[a], b=Flavors[b], c=Flavors[c])
                channel === :zero && continue

                val += f1_x1_vec[a] * f1_x2_vec[b] * D1_vec[c] * H_U[channel]
            end
        end
    end

    αs = αs_func(μ)
    factor = sin(θ)/2 * αs^2/(2*pi*S) * 1/(x1 * x2) * 1/denom

    return factor * val
end

function pp_EEC_U(; pT::Real, θ::Real, η::Real, sqrts::Real, rep::Int)

    S = sqrts^2
    T = -sqrts * pT * exp(-η)
    U = -sqrts * pT * exp( η)

    denom_min = S + U
    denom_min <= 0 && return 0.0
    x2_min = min(max(-T / denom_min, 0.0), 1.0)

    integrand(x2) = pp_EEC_U_integrand(x2=x2, pT=pT, θ=θ, η=η, sqrts=sqrts, rep=rep)

    result = FastGK.integrate(f = integrand, a = x2_min, b = 1.0, rtol = 1e-2, N = 4)
    return result
end

function pp_EEC_U_bin(; pT::Real, θ_lo::Real, θ_hi::Real, η_lo::Real, η_hi::Real, sqrts::Real, rep::Int)

    S = sqrts^2

    x2_lo(v) = begin
        η = v[2]
        T = -sqrts * pT * exp(-η)
        U = -sqrts * pT * exp(η)
        denom = S + U
        denom <= 0 && return 1.0
        return min(max(-T / denom, 0.0), 1.0)
    end
    x2_hi(v) = 1.0

    θ_lo_func(v) = Float64(θ_lo)
    θ_up_func(v) = Float64(θ_hi)

    f_factor(v) = 1.0
    f_NP(v) = 1.0
    f_Pert(v) = pp_EEC_U_integrand(
        x2 = v[1],
        θ = v[2],
        η = v[3],
        pT = pT,
        sqrts = sqrts,
        rep = rep,
    )
    product_rule(pert, NP) = pert * NP

    settings = MultiGKTable.integration_settings(
        f_factor,
        f_Pert,
        f_NP,
        product_rule,
        Float64,
        [(x2_lo, x2_hi), (θ_lo_func, θ_up_func), (Float64(η_lo), Float64(η_hi))],
        [7, 5, 5],
        [5, 3, 3],
        [1e-2, 5e-2, 5e-2],
        [1e-12, 1e-12, 1e-12],
    )

    result, _ = MultiGKTable.integrate(settings)
    return result
end

####################################################################################################
# Asymmetry
####################################################################################################

function pp_EEC_AUT(; pT::Real, θ::Real, η::Real, sqrts::Real, rep::Int)

    num = pp_EEC_T(pT = pT, θ = θ, η = η, sqrts = sqrts, rep = rep)
    den = pp_EEC_U(pT = pT, θ = θ, η = η, sqrts = sqrts, rep = rep)

    return num / den
end

function pp_EEC_AUT_bin(; pT::Real, θ_lo::Real, θ_hi::Real, η_lo::Real, η_hi::Real, sqrts::Real, rep::Int)

    num = pp_EEC_T_bin(pT = pT, θ_lo = θ_lo, θ_hi = θ_hi, η_lo = η_lo, η_hi = η_hi, sqrts = sqrts, rep = rep)
    den = pp_EEC_U_bin(pT = pT, θ_lo = θ_lo, θ_hi = θ_hi, η_lo = η_lo, η_hi = η_hi, sqrts = sqrts, rep = rep)

    return num / den
end

function pp_EEC_AUT_bin_pmap(; pT_array::AbstractVector, θ_lo_array::AbstractVector, θ_hi_array::AbstractVector, η_lo_array::AbstractVector, η_hi_array::AbstractVector, sqrts_array::AbstractVector, rep::AbstractVector)

    n = length(pT_array)
    length(θ_lo_array) == n || throw(ArgumentError("θ_lo_array must have the same length as pT_array"))
    length(θ_hi_array) == n || throw(ArgumentError("θ_hi_array must have the same length as pT_array"))
    length(η_lo_array) == n || throw(ArgumentError("η_lo_array must have the same length as pT_array"))
    length(η_hi_array) == n || throw(ArgumentError("η_hi_array must have the same length as pT_array"))
    length(sqrts_array) == n || throw(ArgumentError("sqrts_array must have the same length as pT_array"))
    length(rep) == n || throw(ArgumentError("rep must have the same length as pT_array"))

    args_vec = collect(zip(pT_array, θ_lo_array, θ_hi_array, η_lo_array, η_hi_array, sqrts_array, Int.(rep)))

    predictions = nothing
    t = @elapsed begin
        predictions = pmap(args_vec; batch_size = 1) do (pT, θ_lo, θ_hi, η_lo, η_hi, sqrts, rep_i)
            pp_EEC_AUT_bin(pT = pT, θ_lo = θ_lo, θ_hi = θ_hi, η_lo = η_lo, η_hi = η_hi, sqrts = sqrts, rep = rep_i)
        end
    end

    return predictions, t
end