function SIDIS_EEC(; chi::Real, x::Real, Q::Real, mu::Real)

    z = (1 - cos(chi)) / 2

    μpdf = Float64(mu)
    μf = max(μpdf, 1.0)
    U = D1_EEC_own_precompute_evolution(Float64(mb), μf)

    D0q = p0 * Q^2 * exp(-z * Q^2 / p1) / (1 + z * Q^2 / p2)
    Dq = U[1, 1] * D0q

    return sin(chi) / 2 * Dq * sum(eq2_vec(μpdf) .* get_f1(x, μpdf))
end

function SIDIS_EEC_batch(; chi_array::AbstractVector, x::Real, Q::Real, mu::Real)
    μpdf = Float64(mu)
    μf = max(μpdf, 1.0)
    U = D1_EEC_own_precompute_evolution(Float64(mb), μf)

    charge_pdf_sum = sum(eq2_vec(μpdf) .* get_f1(x, μpdf))
    Q2 = Float64(Q)^2
    prefactor = charge_pdf_sum * U[1, 1] * p0 * Q2 / 2

    predictions = Vector{Float64}(undef, length(chi_array))
    @inbounds for i in eachindex(chi_array)
        chi = Float64(chi_array[i])
        z = (1 - cos(chi)) / 2
        predictions[i] = sin(chi) * prefactor * exp(-z * Q2 / p1) / (1 + z * Q2 / p2)
    end

    return predictions
end

function SIDIS_EEC_pmap(; chi_array::AbstractVector, x_array::AbstractVector, Q_array::AbstractVector, mu_array::AbstractVector)
    
    n = length(chi_array)
    length(x_array) == n || throw(ArgumentError("x_array must have the same length as chi_array"))
    length(Q_array) == n || throw(ArgumentError("Q_array must have the same length as chi_array"))
    length(mu_array) == n || throw(ArgumentError("mu_array must have the same length as chi_array"))

    args_vec = collect(zip(chi_array, x_array, Q_array, mu_array))

    predictions = nothing
    t = @elapsed begin
        if n > 0 &&
           all(==(x_array[1]), x_array) &&
           all(==(Q_array[1]), Q_array) &&
           all(==(mu_array[1]), mu_array)
            predictions = SIDIS_EEC_batch(
                chi_array = chi_array,
                x = x_array[1],
                Q = Q_array[1],
                mu = mu_array[1],
            )
        else
            predictions = pmap(args_vec; batch_size = 1) do (chi, x, Q, mu)
                SIDIS_EEC(chi = chi, x = x, Q = Q, mu = mu)
            end
        end
    end

    return predictions, t
end
