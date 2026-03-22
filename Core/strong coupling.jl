include("constants.jl")
#include("anomalous dims.jl")

function nf_func(μ)
    if μ > mb
        value = 5
    elseif μ > mc
        value = 4
    else
        value = 3
    end
    return value
end

function m_func(nf)
    if nf == 4
        value = mc
    elseif nf == 5
        value = mb
    end
    return value
end  

function ΛQCD_func(nf)
    if nf == 3
        value = ΛQCD3 
    elseif nf == 4
        value = ΛQCD4
    elseif nf == 5
        value = ΛQCD5
    end
    return value
end 

function αs_func(μ)
    return get_alphas(0, μ)
end

function αs_func_old(μ) 

    order = 2

    nf = nf_func(μ)

    β0 = β0_func(nf)
    β1 = β1_func(nf)
    β2 = β2_func(nf)
    β3 = β3_func(nf)
    β4 = β4_func(nf)

    Λ = ΛQCD_func(nf)   
    x = log(μ^2/Λ^2)

    order0 = 1/(β0*x)
    order1 = - β1/(β0^3) * log(x)/(x^2) 
    order2 = β1^2/(β0^5) * (log(x)^2 - log(x) - 1)/(x^3) + β2/(β0^4)/x^3

     # https://pdg.lbl.gov/2021/reviews/rpp2021-rev-qcd.pdf

    order3 = (β1^3 * (-2*log(x)^3 + 5*log(x)^2 + 4*log(x) - 1) - 6*β0*β1*β2*log(x) + β0^2*β3
                 ) / (2 * β0^7 * x^4)

    order4 = (β1^4 * (6*log(x)^4 - 26*log(x)^3 - 9*log(x)^2 + 24*log(x) + 7) + 18*β0*β1^2*β2 * (2*log(x)^2 - log(x) - 1) + 
                 2*β0^2 * (5*β2 + β0*β4) - β0^2 * β1 * β3 * (12*log(x) + 1)
                 ) / (6 * β0^9 * x^5)

    if order == 0
        value = order0
    elseif order == 1
        value = order0 + order1
    elseif order == 2
        value = order0 + order1 + order2
    elseif order == 3
        value = order0 + order1 + order2 + order3
    elseif order == 4
        value = order0 + order1 + order2 + order3 + order4
    end

    return 4*π*value
end

function alpha_qed(Q)

    mtau  = 1.777

    αem_mc    = 0.007476296864
    αem_mtau  = 0.007496122052
    αem_mb    = 0.007570302837
    αem_mz    = 0.007874015748

    # electric‐charge sums
    eq23 = 6.0/9.0
    eq24 = 10.0/9.0
    eq25 = 11.0/9.0

    if Q <= mc
        nl  = 2
        b0  = (nl + 3*eq23)/(3*π)
        αem = αem_mc / (1 - b0*αem_mc * log(Q^2/mc^2))
    elseif Q <= mtau
        nl  = 2
        b0  = (nl + 3*eq24)/(3*π)
        αem = αem_mc / (1 - b0*αem_mc * log(Q^2/mc^2))
    elseif Q <= mb
        nl  = 3
        b0  = (nl + 3*eq24)/(3*π)
        αem = αem_mtau / (1 - b0*αem_mtau * log(Q^2/mtau^2))
    elseif Q <= mz
        nl  = 3
        b0  = (nl + 3*eq25)/(3*π)
        αem = αem_mb / (1 - b0*αem_mb * log(Q^2/mb^2))
    else
        nl  = 3
        b0  = (nl + 3*eq25)/(3*π)
        αem = αem_mz / (1 - b0*αem_mz * log(Q^2/mz^2))
    end

    return αem
end

@inline function eq2_vec(μ::Real)
    if μ > mb
        return [1/9, 4/9, 1/9, 1/9, 4/9, 0.0, 4/9, 1/9, 1/9, 4/9, 1/9]
    elseif μ > mc
        return [0.0, 4/9, 1/9, 1/9, 4/9, 0.0, 4/9, 1/9, 1/9, 4/9, 0.0]
    else
        return [0.0, 0.0, 1/9, 1/9, 4/9, 0.0, 4/9, 1/9, 1/9, 0.0, 0.0]
    end
end