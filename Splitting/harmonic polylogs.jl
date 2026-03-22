#https://arxiv.org/pdf/hep-ph/9905237
#Note that the paper's expression for H10 is wrong

using PolyLog

function H0_func(x)
    value = log(x)
    return value
end

function H1_func(x) #
    value = -log(1-x)
    return value
end

function H2_func(x)
    value = reli(2,x)
    return value
end

function Hm10_func(x)
    value = log(x)*log(1+x) + reli(2,-x)
    return value
end

function H00_func(x)
    value = 1/2*(log(x))^2
    return value
end

function H10_func(x) #
    value = reli(2,(1-x)) - z2
    return value
end

function H11_func(x) #
    value = 1/2*(log(1-x))^2
    return value
end