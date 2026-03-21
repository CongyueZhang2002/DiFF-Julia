using DifferentialEquations

function get_f1(x::Real, μ::Real)
    pdfs = get_pdfs(0, Float64(x), Float64(μ))
    return [pdfs[i] / x for i in FlavorIndices]
end

const D1_EEC_own_evolution_cache = Dict{Tuple{Float64, Float64}, Matrix{Float64}}()

function D1_EEC_own_precompute_evolution(μi::Real, μf::Real)
    key = (Float64(μi), Float64(μf))

    return get!(D1_EEC_own_evolution_cache, key) do
        function f!(du, u, _, μ)
            nf = nf_func(μ)
            γ0 = [
                γqq0_func(nf) γqg0_func(nf)
                γgq0_func(nf) γgg0_func(nf)
            ]

            du .= (-αs_func(μ) / (2 * π * μ)) .* (transpose(γ0) * u)
        end

        tstops = [t for t in (Float64(mc), Float64(mb)) if min(key[1], key[2]) < t < max(key[1], key[2])]

        evolve(u0) = begin
            problem = ODEProblem(f!, u0, key)
            sol = solve(problem, Tsit5(); tstops = tstops, reltol = 1e-10, abstol = 1e-12)
            copy(sol.u[end])
        end

        return hcat(
            evolve([1.0, 0.0]),
            evolve([0.0, 1.0]),
        )
    end
end

function D1_EEC_own_evolution(D0, μi::Real, μf::Real)
    U = D1_EEC_own_precompute_evolution(μi, μf)
    return U * Float64.(D0)
end

function D1_EEC_own(z::Real, Q::Real, μ::Real)
    
    D0 = [
        p0 * Q^2 * exp(-z * Q^2 / p1) / (1 + z * Q^2 / p2),
        0.0,
    ]

    μf = max(Float64(μ), 1.0)
    μ0 = Float64(mb)
    return D1_EEC_own_evolution(D0, μ0, μf)
end

function D1_EEC_own_vec(z::Real, Q::Real, μ::Real)
    Dq, Dg = D1_EEC_own(z, Q, μ)
    return [
        Dq, Dq, Dq, Dq, Dq, 
        Dg, 
        Dq, Dq, Dq, Dq, Dq
    ]
end
