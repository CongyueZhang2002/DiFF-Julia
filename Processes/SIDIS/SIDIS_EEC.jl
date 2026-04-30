sqrts_EIC = 100.0
const SIDIS_EEC_min_μ = 1.0

SIDIS_EEC_safe_μ(μ::Real) = max(μ, SIDIS_EEC_min_μ)

function SIDIS_EEC_y(; x::Real, Q::Real)
    y = Q^2 / (x * sqrts_EIC^2)
    if !(isfinite(y) && 0.0 < y < 1.0)
        throw(DomainError(y, "SIDIS_EEC requires 0 < y < 1 for y = Q^2 / (x * sqrts_EIC^2); got x=$(x), Q=$(Q), sqrts_EIC=$(sqrts_EIC)"))
    end
    return y
end

function SIDIS_EEC_unnormalized(; chi::Real, x::Real, Q::Real, μ::Real)

    μ = SIDIS_EEC_safe_μ(μ)
    y = SIDIS_EEC_y(x=x, Q=Q)
    z = (1 - cos(chi)) / 2

    Dq, _ = D1_EEC_own(z, Q, μ)

    prefactor = 2*pi*alpha_qed(Q)^2/Q^2
    structure = (1+(1-y)^2)/y * Dq * sum(eq2_vec(μ) .* get_f1(x, μ))

    return sin(chi) / 2 * prefactor * structure
end

function SIDIS_EEC(; chi::Real, Q::Real, μ::Real, x=nothing)

    μ = SIDIS_EEC_safe_μ(μ)
    z = (1 - cos(chi)) / 2

    Dq, _ = D1_EEC_own(z, Q, μ)

    return sin(chi) / 2 * Dq
end

function SIDIS_EEC_batch(; chi_array::AbstractVector, Q::Real, μ::Real, x=nothing)

    μ = SIDIS_EEC_safe_μ(μ)
    U = D1_EEC_own_precompute_evolution(mb, μ)[1, 1]

    predictions = Vector{Float64}(undef, length(chi_array))
    @inbounds for i in eachindex(chi_array)
        chi = chi_array[i]
        z = (1 - cos(chi)) / 2
        predictions[i] = sin(chi)/ 2 * D1_EEC_own_initial_quark(z, Q) * U
    end

    return predictions
end

function SIDIS_EEC_pmap(; chi_array::AbstractVector, Q_array::AbstractVector, μ_array::AbstractVector, x_array=nothing)
    
    n = length(chi_array)
    length(Q_array) == n || throw(ArgumentError("Q_array must have the same length as chi_array"))
    length(μ_array) == n || throw(ArgumentError("μ_array must have the same length as chi_array"))
    if x_array !== nothing
        length(x_array) == n || throw(ArgumentError("x_array must have the same length as chi_array"))
    end

    args_vec = collect(zip(chi_array, Q_array, μ_array))

    predictions = nothing
    t = @elapsed begin
        if n > 0 &&
           all(==(Q_array[1]), Q_array) &&
           all(==(μ_array[1]), μ_array)
            predictions = SIDIS_EEC_batch(
                chi_array = chi_array,
                Q = Q_array[1],
                μ = μ_array[1],
            )
        else
            predictions = pmap(args_vec; batch_size = 1) do (chi, Q, μ)
                SIDIS_EEC(chi = chi, Q = Q, μ = μ)
            end
        end
    end

    return predictions, t
end
