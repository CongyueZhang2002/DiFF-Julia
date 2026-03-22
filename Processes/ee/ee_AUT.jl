function EE_EEC_T(; chi::Real, chibar::Real, Q::Real, mu::Real, rep::Int)

    z = (1 - cos(chi)) / 2
    zbar = (1 - cos(chibar)) / 2

    H1a_vec = H1a_EEC_vec(zchi = z, Q = Q, mu = mu, rep = rep)
    H1a_bar_vec = H1a_EEC_vec(zchi = zbar, Q = Q, mu = mu, rep = rep)

    return sin(chi) / 2 * sin(chibar) / 2 * sum(eq2_vec(mu) .* H1a_vec .* H1a_bar_vec)
end

function EE_EEC_U(; chi::Real, chibar::Real, Q::Real, mu::Real, rep::Int)

    z = (1 - cos(chi)) / 2
    zbar = (1 - cos(chibar)) / 2

    D1_vec = D1_EEC_vec(zchi = z, Q = Q, mu = mu, rep = rep)
    D1_bar_vec = D1_EEC_vec(zchi = zbar, Q = Q, mu = mu, rep = rep)

    return sin(chi) / 2 * sin(chibar) / 2 * sum(eq2_vec(mu) .* D1_vec .* D1_bar_vec)
end

function EE_EEC_AUT(; chi::Real, chibar::Real, Q::Real, mu::Real, rep::Int)

    num = EE_EEC_T(chi = chi, chibar = chibar, Q = Q, mu = mu, rep = rep)
    den = EE_EEC_U(chi = chi, chibar = chibar, Q = Q, mu = mu, rep = rep)

    return num / den
end

function EE_EEC_AUT_pmap(; chi_array::AbstractVector, chibar_array::AbstractVector, Q_array::AbstractVector, mu_array::AbstractVector, rep::AbstractVector)
    
    n = length(chi_array)
    length(chibar_array) == n || throw(ArgumentError("chibar_array must have the same length as chi_array"))
    length(Q_array) == n || throw(ArgumentError("Q_array must have the same length as chi_array"))
    length(mu_array) == n || throw(ArgumentError("mu_array must have the same length as chi_array"))
    length(rep) == n || throw(ArgumentError("rep must have the same length as chi_array"))

    args_vec = collect(zip(chi_array, chibar_array, Q_array, mu_array, Int.(rep)))

    predictions = nothing
    t = @elapsed begin
        predictions = pmap(args_vec; batch_size = 1) do (chi, chibar, Q, mu, rep_i)
            EE_EEC_AUT(chi = chi, chibar = chibar, Q = Q, mu = mu, rep = rep_i)
        end
    end

    return predictions, t
end