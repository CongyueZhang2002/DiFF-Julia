const M1 = 0.13957
const M2 = 0.13957
const NG = 10
const T = Float64[
    -0.97390652851717174,
    -0.86506336668898454,
    -0.67940956829902444,
    -0.43339539412924721,
    -0.14887433898163122,
    0.14887433898163122,
    0.43339539412924721,
    0.67940956829902444,
    0.86506336668898454,
    0.97390652851717174,
]
const W = Float64[
    0.066671344308688069,
    0.14945134915058036,
    0.21908636251598201,
    0.26926671930999652,
    0.29552422471475298,
    0.29552422471475298,
    0.26926671930999652,
    0.21908636251598201,
    0.14945134915058036,
    0.066671344308688069,
]

@inline function get_mu(Q, scale::AbstractString)
    c = 0.5
    if scale == "Q/2"
        return c * Q / 2
    elseif scale == "2Q"
        return c * 2 * Q
    elseif scale == "Q"
        return c * Q
    else
        throw(ArgumentError("scale must be one of \"Q/2\", \"Q\", or \"2Q\""))
    end
end

@inline function get_DiFF(kind::AbstractString, z, Mh, mu, rep)
    values = JAM_DiFF_grid(kind = kind, z = z, Mh = Mh, mu = mu, id = rep)

    if kind == "H1a"
        H1angle_u = values[1]
        return Float64[
            0.0,
            0.0,
            0.0,
            H1angle_u,
            -H1angle_u,
            0.0,
            H1angle_u,
            -H1angle_u,
            0.0,
            0.0,
            0.0,
        ]
    elseif kind == "D1"
        D1_u = values[1]
        D1_s = values[2]
        D1_c = values[3]
        D1_b = values[4]
        D1_g = values[5]
        return Float64[
            D1_b,
            D1_c,
            D1_s,
            D1_u,
            D1_u,
            D1_g,
            D1_u,
            D1_u,
            D1_s,
            D1_c,
            D1_b,
        ]
    else
        throw(ArgumentError("kind must be \"D1\" or \"H1a\""))
    end
end

@inline get_H1a(z, Mh, mu, rep) = get_DiFF("H1a", z, Mh, mu, rep)
@inline get_D1(z, Mh, mu, rep) = get_DiFF("D1", z, Mh, mu, rep)

function get_DiFF_EEC(kind::AbstractString, zchi, Q, scale::AbstractString, rep)
    mu = max(get_mu(Q, scale), 1.0)

    x_min = 0.0
    x_max = 1.0
    y_min = -1.0
    y_max = 1.0

    if kind == "D1"
        a = a_D1
        b = b_D1
    elseif kind == "H1a"
        a = a_H1a
        b = b_H1a
    end

    scriptD = zeros(Float64, 11)
    indices = kind == "H1a" ? (7,) : (6, 7, 9, 10, 11)

    for ip in indices
        integral_sum = 0.0

        for i in 1:NG
            xi = (x_max - x_min) / 2 * T[i] + (x_max + x_min) / 2

            inner_sum = 0.0
            for j in 1:NG
                yj = (y_max - y_min) / 2 * T[j] + (y_max + y_min) / 2

                term1 = Q^2 * xi^2 * (1 - yj^2) * zchi / 4
                term2 = 2 * (M1^2 / (1 + yj) + M2^2 / (1 - yj))
                Mh = sqrt(term1 + term2)

                diff_value = get_DiFF(kind, xi, Mh, mu, rep)[ip]
                diff_xi_zeta_Mh =
                    0.5 * diff_value +
                    a * diff_value * yj +
                    b * diff_value * 0.5 * (3 * yj^2 - 1)

                inner_sum += W[j] * xi^4 * (1 - yj^2)^2 / Mh * diff_xi_zeta_Mh
            end

            integral_sum += W[i] * inner_sum * ((y_max - y_min) / 2)
        end

        scriptD[ip] = Q^2 / 32 * integral_sum * ((x_max - x_min) / 2)
    end

    if kind == "H1a"
        scriptH_u = scriptD[7]
        return Float64[
            0.0,
            0.0,
            0.0,
            scriptH_u,
            -scriptH_u,
            0.0,
            scriptH_u,
            -scriptH_u,
            0.0,
            0.0,
            0.0,
        ]
    else
        return Float64[
            scriptD[11],
            scriptD[10],
            scriptD[9],
            scriptD[7],
            scriptD[7],
            scriptD[6],
            scriptD[7],
            scriptD[7],
            scriptD[9],
            scriptD[10],
            scriptD[11],
        ]
    end
end

@inline get_H1a_EEC(zchi, Q, scale::AbstractString, rep) =
    get_DiFF_EEC("H1a", zchi, Q, scale, rep)

@inline get_D1_EEC(zchi, Q, scale::AbstractString, rep) =
    get_DiFF_EEC("D1", zchi, Q, scale, rep)
