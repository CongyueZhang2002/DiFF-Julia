function SIDIS_EEC(; chi::Real, x::Real, Q::Real, μ::Real)

    z = (1 - cos(chi)) / 2

    Dq, _ = D1_EEC_own(z, Q, μ)

    return sin(chi) / 2 * Dq * sum(eq2_vec(μ) .* get_f1(x, μ))
end

function SIDIS_EEC_batch(; chi_array::AbstractVector, x::Real, Q::Real, μ::Real)

    U = D1_EEC_own_precompute_evolution(mb, max(μ, 1.0))[1, 1]

    factor = sum(eq2_vec(μ) .* get_f1(x, μ)) * U

    predictions = Vector{Float64}(undef, length(chi_array))
    @inbounds for i in eachindex(chi_array)
        chi = chi_array[i]
        z = (1 - cos(chi)) / 2
        predictions[i] = sin(chi)/ 2 * factor * D1_EEC_own_initial_quark(z, Q)
    end

    return predictions
end

function SIDIS_EEC_pmap(; chi_array::AbstractVector, x_array::AbstractVector, Q_array::AbstractVector, μ_array::AbstractVector)
    
    n = length(chi_array)
    length(x_array) == n || throw(ArgumentError("x_array must have the same length as chi_array"))
    length(Q_array) == n || throw(ArgumentError("Q_array must have the same length as chi_array"))
    length(μ_array) == n || throw(ArgumentError("μ_array must have the same length as chi_array"))

    args_vec = collect(zip(chi_array, x_array, Q_array, μ_array))

    predictions = nothing
    t = @elapsed begin
        if n > 0 &&
           all(==(x_array[1]), x_array) &&
           all(==(Q_array[1]), Q_array) &&
           all(==(μ_array[1]), μ_array)
            predictions = SIDIS_EEC_batch(
                chi_array = chi_array,
                x = x_array[1],
                Q = Q_array[1],
                μ = μ_array[1],
            )
        else
            predictions = pmap(args_vec; batch_size = 1) do (chi, x, Q, μ)
                SIDIS_EEC(chi = chi, x = x, Q = Q, μ = μ)
            end
        end
    end

    return predictions, t
end
