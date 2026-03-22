function EE_EEC(; chi::Real, Q::Real, mu::Real, rep::Int)

    z = (1 - cos(chi)) / 2

    D1_vec = D1_EEC_vec(zchi = z, Q = Q, mu = mu, rep = rep)
    eq2 = eq2_vec(mu)
    charge_norm = sum(eq2)

    Dq = sum(eq2 .* D1_vec) / charge_norm

    return sin(chi) / 4 * Dq
end

function EE_EEC_pmap(; chi_array::AbstractVector, Q_array::AbstractVector, mu_array::AbstractVector, rep::AbstractVector)
    
    n = length(chi_array)
    length(Q_array) == n || throw(ArgumentError("Q_array must have the same length as chi_array"))
    length(mu_array) == n || throw(ArgumentError("mu_array must have the same length as chi_array"))
    length(rep) == n || throw(ArgumentError("rep must have the same length as chi_array"))

    args_vec = collect(zip(chi_array, Q_array, mu_array, Int.(rep)))

    predictions = nothing
    t = @elapsed begin
        predictions = pmap(args_vec; batch_size = 1) do (chi, Q, mu, rep_i)
            EE_EEC(chi = chi, Q = Q, mu = mu, rep = rep_i)
        end
    end

    return predictions, t
end
