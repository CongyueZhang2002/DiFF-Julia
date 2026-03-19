# Precomputed anomalous dims

# 11/3*CA - 4/3*TF*nf
function β0_func(nf::Integer)
    if nf == 3
        return 9.0
    elseif nf == 4
        return 8.333333333333334
    elseif nf == 5
        return 7.666666666666667
    else
        error("β0_func: unsupported nf = $nf")
    end
end

# 34/3*CA^2 - 20/3*CA*TF*nf - 4*CF*TF*nf
function β1_func(nf::Integer)
    if nf == 3
        return 64.0
    elseif nf == 4
        return 51.333333333333336
    elseif nf == 5
        return 38.66666666666667
    else
        error("β1_func: unsupported nf = $nf")
    end
end

#((2857/54)*CA^3 + (2*CF^2 - (205/9)*CF*CA - (1415/27)*CA^2)*TF*nf + ((44/9)*CF + (158/27)*CA)*TF^2*nf^2)
function β2_func(nf::Integer)
    if nf == 3
        return 643.8333333333334
    elseif nf == 4
        return 406.35185185185196
    elseif nf == 5
        return 180.90740740740756
    else
        error("β2_func: unsupported nf = $nf")
    end
end

# https://arxiv.org/pdf/hep-ph/9701390
function β3_func(nf)
    if nf == 3
        return 1.209037813024e+04
    elseif nf == 4
        return 8.035186419171e+03
    elseif nf == 5
        return 4.826156328096e+03
    else
        error("β3_func: unsupported nf = $nf")
    end
end

# https://arxiv.org/pdf/1606.08659 (uses alphs/pi convention so have to renomralize)
function β4_func(nf)
    if nf == 3
        return 1.303779068020e+05
    elseif nf == 4
        return 5.831055395044e+04
    elseif nf == 5
        return 1.547061222594e+04
    else
        error("β4_func: unsupported nf = $nf")
    end
end

# 4*CF
function Γ0_func(nf::Integer)
    if nf == 3
        return 5.333333333333333
    elseif nf == 4
        return 5.333333333333333
    elseif nf == 5
        return 5.333333333333333
    else
        error("Γ0_func: unsupported nf = $nf")
    end
end

# https://arxiv.org/pdf/1911.10174

# 4*CF*((67/9 - π^2/3)*CA - 20/9*TF*nf)
function Γ1_func(nf::Integer)
    if nf == 3
        return 48.69544319419008
    elseif nf == 4
        return 42.76951726826415
    elseif nf == 5
        return 36.84359134233823
    else
        error("Γ1_func: unsupported nf = $nf")
    end
end

# 4 * CF * (CA^2 * (245/6 - (134*π^2)/27 + (11*π^4)/45 + (22/3) * z3) + CA * TF * nf * ( -  418/27  + (40*π^2)/27  - (56/3) * z3) + CF * TF * nf * (  -  55/3  + 16 * z3) - 16/27 * TF^2 * nf^2)
function Γ2_func(nf::Integer)
    if nf == 3
        return 618.2248693799629
    elseif nf == 4
        return 429.5065747550467
    elseif nf == 5
        return 239.20803321655015
    else
        error("Γ2_func: unsupported nf = $nf")
    end
end

function Γ3_func(nf)
    if nf == 3
        return 7.035152974055e+03
    elseif nf == 4
        return 3.353354170411e+03
    elseif nf == 5
        return 1.412460850831e+02
    else
        error("Γ3_func: unsupported nf = $nf")
    end
end

function Γ4_func(nf::Integer)
    if nf == 3
        return 0.21
    elseif nf == 4
        return 0.21
    elseif nf == 5
        return 0.21
    else
        error("Γ4_func: unsupported nf = $nf")
    end
end

function Γ_func(; μ::Float64, order::Int64)

    as = αs_func(μ)/(4π)
    nf = nf_func(μ)

    Γ0 = Γ0_func(nf)
    Γ1 = Γ1_func(nf)
    Γ2 = Γ2_func(nf)
    Γ3 = Γ3_func(nf)
    Γ4 = Γ4_func(nf)

    if order == 1
        value = as*Γ0 + as^2*Γ1
    elseif order == 2
        value = as*Γ0 + as^2*Γ1 + as^3*Γ2
    elseif order == 3
        value = as*Γ0 + as^2*Γ1 + as^3*Γ2 + as^4*Γ3
    elseif order == 4
        value = as*Γ0 + as^2*Γ1 + as^3*Γ2 + as^4*Γ3 + as^5*Γ4
    else
        error("Γ_func: unsupported order = $order")
    end
end