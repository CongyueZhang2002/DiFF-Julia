using QCDNUM
#include("test.jl")

function FF_func(; z::Float64, μ::Float64, flavor::String)
    itype = 3 # FF

    if μ^2 < 1.01
        μ_safe = sqrt(1.01)
    else
        μ_safe = μ 
    end

    FFs = QCDNUM.allfxq(itype, z, μ_safe^2, 0, 1)

    if flavor == "bb"
        i = 12
    elseif flavor == "cb"
        i = 11
    elseif flavor == "sb"
        i = 10
    elseif flavor == "ub"
        i = 8
    elseif flavor == "db"
        i = 9
    elseif flavor == "g"
        i = 7
    elseif flavor == "d"
        i = 5
    elseif flavor == "u"
        i = 6
    elseif flavor == "s"
        i = 4
    elseif flavor == "c"
        i = 3
    elseif flavor == "b"
        i = 2
    end

    value = FFs[i]/z

    return value
end

#display(FF_func(z=0.001, μ=2.0, flavor="u"))