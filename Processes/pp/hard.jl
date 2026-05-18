# Unpolarized hard functions for 2->2 from RevModPhys.59.465.

@inline H_U_qqp_qqp_func(; s::Float64, t::Float64, u::Float64) = (
    4/9 * (s^2 + u^2) / t^2
)

@inline H_U_qq_qq_func(; s::Float64, t::Float64, u::Float64) = (
    4/9 * ((s^2 + u^2) / t^2 + (s^2 + t^2) / u^2) -
    8/27 * s^2 / (t * u)
)

@inline H_U_qqb_qpqpb_func(; s::Float64, t::Float64, u::Float64) = (
    4/9 * (t^2 + u^2) / s^2
)

@inline H_U_qqb_qqb_func(; s::Float64, t::Float64, u::Float64) = (
    4/9 * ((s^2 + u^2) / t^2 + (t^2 + u^2) / s^2) -
    8/27 * u^2 / (s * t)
)

@inline H_U_gq_gq_func(; s::Float64, t::Float64, u::Float64) = (
    (s^2 + u^2) / t^2 -
    4/9 * (s^2 + u^2) / (s * u)
)

@inline H_U_qqb_gg_func(; s::Float64, t::Float64, u::Float64) = (
    32/27 * (t/u + u/t) -
    8/3 * (t^2 + u^2) / s^2
)

@inline H_U_gg_qqb_func(; s::Float64, t::Float64, u::Float64) = (
    1/6 * (t/u + u/t) -
    3/8 * (t^2 + u^2) / s^2
)

@inline H_U_gg_gg_func(; s::Float64, t::Float64, u::Float64) = (
    9/2 * (3 - (t*u)/s^2 - (s*u)/t^2 - (s*t)/u^2)
)

@inline H_U_func(; s::Float64, t::Float64, u::Float64) = Dict(

    "qqp_qqp" => H_U_qqp_qqp_func(s=s, t=t, u=u),
    "qq_qq" => H_U_qq_qq_func(s=s, t=t, u=u),
    "qqb_qpqpb" => H_U_qqb_qpqpb_func(s=s, t=t, u=u),
    "qqb_qqb" => H_U_qqb_qqb_func(s=s, t=t, u=u),
    "gq_gq" => H_U_gq_gq_func(s=s, t=t, u=u),
    "qqb_gg" => H_U_qqb_gg_func(s=s, t=t, u=u),
    "gg_qqb" => H_U_gg_qqb_func(s=s, t=t, u=u),
    "gg_gg" => H_U_gg_gg_func(s=s, t=t, u=u),

    "gq_gq_cross" => H_U_gq_gq_func(s=s, t=u, u=t),
    "qqp_qqp_cross" => H_U_qqp_qqp_func(s=s, t=u, u=t),
    "qqb_qqb_cross" => H_U_qqb_qqb_func(s=s, t=u, u=t),

    "zero" => 0.0
)

@inline qb_func(flavor::String) = (
    if flavor == "bb"
        "b"
    elseif flavor == "cb"
        "c"
    elseif flavor == "sb"
        "s"
    elseif flavor == "db"
        "d"
    elseif flavor == "ub"
        "u"
    elseif flavor == "g"
        "g"
    elseif flavor == "u"
        "ub"
    elseif flavor == "d"
        "db"
    elseif flavor == "s"
        "sb"
    elseif flavor == "c"
        "cb"
    elseif flavor == "b"
        "bb"
    else
        throw(ArgumentError("unknown flavor: $flavor"))
    end
)

@inline if_qp(flavor1::String, flavor2::String) = (
    flavor2 != flavor1 && flavor2 != "g" && flavor2 != qb_func(flavor1)
)

@inline function H_U_channel(; a::String, b::String, c::String)
    # a+b->c+d
    # c is the fragmenting parton.

    if a == "g"
        if b == "g"
            if c == "g"
                return "gg_gg"
            elseif c != "g"
                return "gg_qqb"
            else
                return "zero"
            end
        elseif b != "g"
            if c == "g"
                return "gq_gq"
            elseif c == b
                return "gq_gq_cross"   
            else
                return "zero"
            end
        end
    else
        if b == "g"
            if c == a
                return "gq_gq"
            elseif c == "g"
                return "gq_gq_cross"   
            else
                return "zero"
            end
        elseif b == a
            if c == a
                return "qq_qq"
            else
                return "zero"
            end
        elseif if_qp(a, b)
            if c == a
                return "qqp_qqp"
            elseif c == b
                return "qqp_qqp_cross" 
            else
                return "zero"
            end
        elseif b == qb_func(a)
            if c == a
                return "qqb_qqb"
            elseif c == b
                return "qqb_qqb_cross" 
            elseif c == "g"
                return "qqb_gg"
            elseif if_qp(a, c)
                return "qqb_qpqpb"
            else
                return "zero"
            end
        else
            return "zero"
        end
    end
end

# Transversely polarized hard functions for 2->2 from 0709.3272 and AN_theory0_hadinjet.py

@inline H_T_qg_qg_func(; s::Float64, t::Float64, u::Float64) = (
    8/9 - 2 * s * u / t^2
)

@inline H_T_qq_qq_func(; s::Float64, t::Float64, u::Float64) = (
    8/9 * (-s * u / t^2 + 1/3 * s / t)
)

@inline H_T_qqb_qqb_func(; s::Float64, t::Float64, u::Float64) = (
    8/9 * (-s * u / t^2 + 1/3 * u / t)
)

@inline H_T_qqb_qbq_func(; s::Float64, t::Float64, u::Float64) = (
    -8/27
)

@inline H_T_qqp_qqp_func(; s::Float64, t::Float64, u::Float64) = (
    -8/9 * s * u / t^2
)

@inline H_T_func(; s::Float64, t::Float64, u::Float64) = Dict(
    "qg_qg" => H_T_qg_qg_func(s=s, t=t, u=u),
    "qq_qq" => H_T_qq_qq_func(s=s, t=t, u=u),
    "qqb_qqb" => H_T_qqb_qqb_func(s=s, t=t, u=u),
    "qqb_qbq" => H_T_qqb_qbq_func(s=s, t=t, u=u),
    "qqp_qqp" => H_T_qqp_qqp_func(s=s, t=t, u=u),
    "zero" => 0.0
)

@inline function H_T_channel(; a::String, b::String, c::String)
    # a+b->c+d
    # a is the polarized-side parton.
    # c is the fragmenting parton.

    if a == "g"
        return "zero"
    elseif b == a && c == a
        return "qq_qq"
    elseif if_qp(a, b) && c == a
        return "qqp_qqp"
    elseif b == qb_func(a) && c == a
        return "qqb_qqb"
    elseif b == qb_func(a) && c == b
        return "qqb_qbq"
    elseif b == "g" && c == a
        return "qg_qg"
    else
        return "zero"
    end
end