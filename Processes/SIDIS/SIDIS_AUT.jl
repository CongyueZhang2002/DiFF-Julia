function SIDIS_EEC_T(; chi::Real, x::Real, Q::Real, mu::Real, rep::Int)

    z = (1 - cos(chi)) / 2

    H1a_vec = H1a_EEC_vec(zchi = z, Q = Q, mu = mu, rep = rep)
    h1_vec = get_h1(x, mu, rep)

    return - sin(chi) / 2 * sum(eq2_vec(mu) .* H1a_vec .* h1_vec)
end

function SIDIS_EEC_U(; chi::Real, x::Real, Q::Real, mu::Real, rep::Int)

    z = (1 - cos(chi)) / 2

    D1_vec = D1_EEC_vec(zchi = z, Q = Q, mu = mu, rep = rep)
    f1_vec = get_f1(x, mu)

    return sin(chi) / 2 * sum(eq2_vec(mu) .* D1_vec .* f1_vec)
end

function SIDIS_EEC_AUT(; chi::Real, x::Real, Q::Real, mu::Real, rep::Int)

    num = SIDIS_EEC_T(chi=chi, x=x, Q=Q, mu=mu, rep=rep)
    den = SIDIS_EEC_U(chi=chi, x=x, Q=Q, mu=mu, rep=rep)

    return num/den
end

function SIDIS_EEC_AUT_pmap(; chi_array::AbstractVector, x_array::AbstractVector, Q_array::AbstractVector, mu_array::AbstractVector, rep::AbstractVector)
    
    n = length(chi_array)
    length(x_array) == n || throw(ArgumentError("x_array must have the same length as chi_array"))
    length(Q_array) == n || throw(ArgumentError("Q_array must have the same length as chi_array"))
    length(mu_array) == n || throw(ArgumentError("mu_array must have the same length as chi_array"))
    length(rep) == n || throw(ArgumentError("rep must have the same length as chi_array"))

    args_vec = collect(zip(chi_array, x_array, Q_array, mu_array, Int.(rep)))

    predictions = nothing
    t = @elapsed begin
        predictions = pmap(args_vec; batch_size = 1) do (chi, x, Q, mu, rep_i)
            SIDIS_EEC_AUT(chi = chi, x = x, Q = Q, mu = mu, rep = rep_i)
        end
    end

    return predictions, t
end