using DifferentialEquations

function get_f1(x::Real, μ::Real)
    pdfs = get_pdfs(0, Float64(x), Float64(μ))
    return [pdfs[i] for i in FlavorIndices]
end

const D1_EEC_own_evolution_cache = Dict{Tuple{Float64, Float64}, Matrix{Float64}}()

@inline function D1_EEC_own_anomalous_dimension(μ::Real)
    nf = nf_func(μ)
    return [
        γqq0_func(nf) γqg0_func(nf)
        γgq0_func(nf) γgg0_func(nf)
    ]
end

@inline function D1_EEC_own_tstops(μi::Float64, μf::Float64)
    μmin, μmax = minmax(μi, μf)
    return [
        threshold for threshold in (Float64(mc), Float64(mb))
        if μmin < threshold < μmax
    ]
end

function D1_EEC_own_precompute_evolution(μi::Real, μf::Real)
    key = (Float64(μi), Float64(μf))

    return get!(D1_EEC_own_evolution_cache, key) do
        function f!(du, u, _, μ)
            γ0 = D1_EEC_own_anomalous_dimension(μ)
            du .= (-αs_func(μ) / (2 * π * μ)) .* (transpose(γ0) * u)
        end

        problem = ODEProblem(f!, [1.0 0.0; 0.0 1.0], key)
        sol = solve(
            problem,
            Tsit5();
            tstops = D1_EEC_own_tstops(key...),
            reltol = 1e-10,
            abstol = 1e-12,
        )
        return Matrix(sol.u[end])
    end
end

function D1_EEC_own_evolution(D0, μi::Real, μf::Real)
    U = D1_EEC_own_precompute_evolution(μi, μf)
    return U * Float64.(D0)
end

@inline function D1_EEC_own_initial_quark(z::Real, Q::Real)
    Q2 = Float64(Q)^2
    zQ2 = Float64(z) * Q2
    return p0 * Q2 * exp(-zQ2 / p1) / (1 + zQ2 / p2)
end

function D1_EEC_own(z::Real, Q::Real, μ::Real)
    μf = max(Float64(μ), 1.0)
    D0 = [D1_EEC_own_initial_quark(z, Q), 0.0]
    return D1_EEC_own_evolution(D0, Float64(mb), μf)
end

function D1_EEC_own_vec(z::Real, Q::Real, μ::Real)
    Dq, Dg = D1_EEC_own(z, Q, μ)
    return [fill(Dq, 5); Dg; fill(Dq, 5)]
end
