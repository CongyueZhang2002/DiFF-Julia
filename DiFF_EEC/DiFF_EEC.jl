const Flavors = ("bb", "cb", "sb", "db", "ub", "g", "u", "d", "s", "c", "b")
const FlavorIndices = (-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)

function get_f1(x::Real, μ::Real)
    pdfs = get_pdfs(0, Float64(x), Float64(μ))
    return [pdfs[i] / x for i in FlavorIndices]
end

function get_h1(x::Real, μ::Real, rep::Int)
    set_lhapdf(1, pdf_dict_array[2]["pdfset_name"], rep)
    pdfs = get_pdfs(1, Float64(x), Float64(μ))
    return [pdfs[i] / x for i in FlavorIndices]
end

@inline function DiFF_EEC(; kind::String, zchi::Real, Q::Real, mu::Real, rep::Int, nth::Int)

    if kind == "D1"
        NP_a, NP_b = a_D1, b_D1
    elseif kind == "H1a"
        NP_a, NP_b = a_H1a, b_H1a
    else
        throw(ArgumentError("kind must be \"D1\" or \"H1a\""))
    end

    integrand_2(x,y) = begin

        Mh = sqrt(
              Q^2 * x^2 * (1 - y^2) * zchi/4 
            + 2 * mpi^2 * (1/(1 + y) + 1/(1 - y))
        )

        JAM_DiFF = JAM_DiFF_grid(kind = kind, z = x, Mh = Mh, mu = mu, rep = rep)[nth]

        value = JAM_DiFF * x^4 * (1 - y^2)^2 / Mh * (
            0.5 + NP_a * y + NP_b * 0.5 * (3 * y^2 - 1)
        )

        return value
    end

    integrand_1(x) = FastGK.integrate(f = y -> integrand_2(x, y), a = -1.0, b = 1.0, rtol = 1e-2, N = 4)

    value = Q^2 / 32 * FastGK.integrate(f = integrand_1, a = 0.0, b = 1.0, rtol = 1e-2, N = 4)

    return value
end

function DiFF_EEC_vec(; kind::String, zchi::Real, Q::Real, mu::Real, rep::Int)

    if kind == "D1"
        D1_u = DiFF_EEC(kind=kind, zchi=zchi, Q=Q, mu=mu, rep=rep, nth=1)
        D1_s = DiFF_EEC(kind=kind, zchi=zchi, Q=Q, mu=mu, rep=rep, nth=2)
        D1_c = DiFF_EEC(kind=kind, zchi=zchi, Q=Q, mu=mu, rep=rep, nth=3)
        D1_b = DiFF_EEC(kind=kind, zchi=zchi, Q=Q, mu=mu, rep=rep, nth=4)
        D1_g = DiFF_EEC(kind=kind, zchi=zchi, Q=Q, mu=mu, rep=rep, nth=5)
        return [
            D1_b, D1_c, D1_s, D1_u, D1_u,
            D1_g,
            D1_u, D1_u, D1_s, D1_c, D1_b
        ]
    else
        H1a_u = DiFF_EEC(kind=kind, zchi=zchi, Q=Q, mu=mu, rep=rep, nth=1)
        return [
            0.0, 0.0, 0.0, H1a_u, -H1a_u,
            0.0,
            H1a_u,-H1a_u, 0.0, 0.0, 0.0
        ]
    end
end

@inline D1_EEC_vec(; zchi::Real, Q::Real, mu::Real, rep::Int) =
    DiFF_EEC_vec(kind="D1", zchi=zchi, Q=Q, mu=mu, rep=rep)

@inline H1a_EEC_vec(; zchi::Real, Q::Real, mu::Real, rep::Int) =
    DiFF_EEC_vec(kind="H1a", zchi=zchi, Q=Q, mu=mu, rep=rep)