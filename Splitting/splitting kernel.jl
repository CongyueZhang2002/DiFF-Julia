#https://www.sciencedirect.com/science/article/pii/S0550321306004330
#timelike singlet splitting functions

using PolyLog
include("harmonic polylogs.jl")
include("../Core/constants.jl")

#Pqq-----------------------------------------------------------------------------------------------

function Pqq0_Reg_func(x,nf)

    value =  - (8*(1 + x))/3

    return value
end

function Pqq0_Delta_func(x,nf)

    value = 4

    return value
end

function Pqq0_D0_func(x,nf)

    value = 16/3

    return value
end

function Pqq1_Reg_func(x,nf) #Pns+ + Pps

    H0 = H0_func(x)
    H2 = H2_func(x)
    Hm10 = Hm10_func(x)
    H00 = H00_func(x)   
    H10 = H10_func(x)

    value = (8*(27 - 146*nf + 8*π^2 - (40*nf)/x - 429*x + 94*nf*x + 10*π^2*x + 112*nf*x^2 + (2*π^2)/(1 + x) - 3*(57 + 41*x + 4*nf*(7 + 13*x + 4*x^2))*H0 + 12*(6 + 5*x - (1 + x)^(-1) + 3*nf*(1 + x))*H00 + 
    48*H10 + 48*x*H10 + 48*H2 + 48*x*H2 - 12*Hm10 + 12*x*Hm10 + (24*Hm10)/(1 + x)))/27
    
    return value
end

function Pqq1_Delta_func(x,nf)

    value = (-2*(-189 - 84*π^2 + nf*(6 + 8*π^2) + 72*z3))/27

    return value
end

function Pqq1_D0_func(x,nf)

    H0 = H0_func(x)
    H2 = H2_func(x)
    H00 = H00_func(x)   
    H10 = H10_func(x)

    value = (-16*(-201 + 10*nf + 9*π^2 + 3*(-45 + 2*nf)*H0 + 42*H00 + 48*H10 + 48*H2))/27

    return value
end

function Pqq2_Reg_func(x,nf)

    L0 = log(x)
    L1 = log(1-x)

    value = (x*(1658.69 + x*(-4249.4 + x*(-1075.3 + 593.9*x))) + nf*(479.869 + x*(-784.620 + x*(-276.030 + x*(1398.94 + x*(-461.914 + 61.284*x))))) + 
        nf^2*(-1.58024 + x*(94.9061 + x*(-150.223 + x*(15.8609 + x*(67.455 + x*(-33.621 + 7.9934*x)))))) + 
        (-39.1111*x + nf*(-28.4444 + x*(78.1368 - 47.322*x + nf*(-1.9000 + x*(-5.96449 + 7.8645*x)))))*L0^3 + 
        (1.58024 + nf*(9.072 + nf*(0.019122 - 0.019122*x) - 9.072*x))*x*L0^4 + 
        L0^2*(-189.37*x + nf*(-14.2222 + x*(437.189 + 9.76079*nf - 425.14*x - 9.76079*nf*x)) + x*(-519.37 + 28.551*nf - 23.102*nf*x)*L1) + 
        L0*(1327.5*x + nf*(324.069 + x*(163.529 - 695.602*x + nf*(62.0586 + (-35.7646 - 26.2939*x)*x))) + 
        x*(-56.9069 - 559.1*x + nf*(45.413 - 2.1031*nf - 10.4510*x + 2.1031*nf*x))*L1) + 
        x*L1*(-707.67 + nf*(54.5598 + nf*(16.611 - 16.611*x) + 8.65*x) + nf*L1*(-9.751 + 1.77799*nf + 9.751*x - 1.77799*nf*x + (-5.926 + 5.926*x)*L1)))/x
    
    return value
end

function Pqq2_Delta_func(x,nf)

    value = 1295.625 + nf*(-173.934 + 1.13066*nf)

    return value 
end

function Pqq2_D0_func(x,nf)

    L0 = log(x)

    value = 1174.89 + (-183.187 - 0.790123*nf)*nf + nf^2*x*L0*(3.95061 + 1.18518*L0)

    return value
end

function Pqq_integrand(; x::Float64, α::Float64, order::Int64, nf::Int64, f::Function)

    Pqq_Reg = Pqq0_Reg_func(x,nf)
    Pqq_Delta = Pqq0_Delta_func(x,nf)
    Pqq_D0_x = Pqq0_D0_func(x,nf)
    Pqq_D0_1 = Pqq0_D0_func(1,nf)
    integrand0 = (
          Pqq_Reg*f(x)
        + Pqq_Delta*f(1) #This works because integration interval is 1
        + (Pqq_D0_x*f(x) - Pqq_D0_1*f(1))/(1-x)
    )

    Pqq_Reg = Pqq1_Reg_func(x,nf)
    Pqq_Delta = Pqq1_Delta_func(x,nf)
    Pqq_D0_x = Pqq1_D0_func(x,nf)
    Pqq_D0_1 = Pqq1_D0_func(1,nf)
    integrand1 = (
          Pqq_Reg*f(x)
        + Pqq_Delta*f(1) #This works because integration interval is 1
        + (Pqq_D0_x*f(x) - Pqq_D0_1*f(1))/(1-x)
    )

    Pqq_Reg = Pqq2_Reg_func(x,nf)
    Pqq_Delta = Pqq2_Delta_func(x,nf)
    Pqq_D0_x = Pqq2_D0_func(x,nf)
    Pqq_D0_1 = Pqq2_D0_func(1,nf)
    integrand2 = (
          Pqq_Reg*f(x)
        + Pqq_Delta*f(1) #This works because integration interval is 1
        + (Pqq_D0_x*f(x) - Pqq_D0_1*f(1))/(1-x)
    )

    #if order == 0
    #    total = (α/(4*π))*integrand0
    #elseif order == 1
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1
    #elseif order == 2
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1 + (α/(4*π))^3*integrand2
    #end

    if order == 0
        total = integrand0
    elseif order == 1
        total = integrand1
    elseif order == 2
        total = integrand2
    end

    return total
end

#Pqg-----------------------------------------------------------------------------------------------

function Pqg0_func(x,nf)

    value = nf*(2 - 4x + 4x^2)

    return value
end

function Pqg1_func(x,nf)

    H0 = H0_func(x)
    H1 = H1_func(x)
    H2 = H2_func(x)
    Hm10 = Hm10_func(x)
    H00 = H00_func(x)   
    H10 = H10_func(x)
    H11 = H11_func(x)

    value = (-4*nf*(60 + x*(33 + 3*(49 - 138*x)*x + 2*nf*(5 + 4*(-1 + x)*x) + 2*π^2*(2 + x*(5 + 4*x))) + 6*x*(11 + nf + 2*nf*(-1 + x)*x + x*(47 + 2*x))*H0 + 6*x*(-7 - 58*x + 8*x^2)*H00 + 
    3*x*((21 + 22*(-1 + x)*x + nf*(-2 - 4*(-1 + x)*x))*H1 - 2*(1 + 2*(-1 + x)*x)*(6*H10 + 5*H11 - 14*H2) + 18*(1 + 2*x*(1 + x))*Hm10)))/(9*x)

    return value
end

function Pqg2_func(x,nf)
    
    L0 = log(x)
    L1 = log(1-x)

    value = (
        nf * (
            100/27 * L1^4 + 350/9 * L1^3 + 263.07 * L1^2 + 693.84 * L1 + 603.71 -
            882.48 * x + 4723.2 * x^2 - 4745.8 * x^3 - 175.28 * x^4 -
            L0 * L1 * (1809.4 + 107.59*x) - 885.5 * x * L0^4 +
            1864 * L0 + 1512 * L0^2 + 361.28 * L0^3 + 42.328 * L0^4 +
            1141.7 / x + 675.83 / x * L0 - 64 / x * (L0^2 + L0^3)
        )
        +
        nf^2 * (
            -100/27 * L1^3 - 35.446 * L1^2 - 103.609 * L1 - 113.81 +
            341.26 * x - 853.35 * x^2 + 492.1 * x^3 + 14.803 * x^4 +
            L0 * L1 * (966.96 - 1.593 * L1 - 709.1 * x) -
            333.8 * x * L0^3 + 619.75 * L0 + 255.62 * L0^2 + 21.569 * L0^3 -
            2.8986 / x - 3.1752 / x * L0 - 32 / 27 / x * L0^2
        )
        +
        nf^3 * (
            4 + 6 * (L0 + L1) + (1 - 2 * x + 2 * x^2) *
            (3.8696 + 4 * (L0 + L1) + 3 * (L0 + L1)^2)
        ) * 4 / 9
    )

    return value
end

function Pqg_integrand(; x::Float64, α::Float64, order::Int64, nf::Int64, f::Function)

    Pqg = Pqg0_func(x,nf)
    integrand0 = (
        Pqg*f(x)
    )

    Pqg = Pqg1_func(x,nf)
    integrand1 = (
        Pqg*f(x)
    )

    Pqg = Pqg2_func(x,nf)
    integrand2 = (
        Pqg*f(x)
    )

    #if order == 0
    #    total = (α/(4*π))*integrand0
    #elseif order == 1
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1
    #elseif order == 2
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1 + (α/(4*π))^3*integrand2
    #end

    if order == 0
        total = integrand0
    elseif order == 1
        total = integrand1
    elseif order == 2
        total = integrand2
    end

    return total
end

#Pgq-----------------------------------------------------------------------------------------------

function Pgq0_func(x,nf)

    value = 4/3*(-4 + 4/x + 2x)

    return value
end

function Pgq1_func(x,nf)

    H0 = H0_func(x)
    H1 = H1_func(x)
    H2 = H2_func(x)
    Hm10 = Hm10_func(x)
    H00 = H00_func(x)   
    H10 = H10_func(x)
    H11 = H11_func(x)

    value = (16*(17 + x*(43 + 6*π^2 + (9 - 44*x)*x) + (-54 + x*(40 + x*(83 + 24*x)))*H0 - 2*(36 + x*(14 + 29*x))*H00 - 76*H10 - 20*H11 + 4*H2 + 36*Hm10 + 
    2*x*(5*x*H1 - (-2 + x)*(19*H10 + 5*H11 - H2) + 9*(2 + x)*Hm10)))/(9*x)

    return value
end

function Pgq2_func(x,nf)
    
    L0 = log(x)
    L1 = log(1-x)

    value = (
        400/81*L1^4 + 520/27*L1^3 − 220.13*L1^2 − 152.6*L1 + 272.93 − 7188.7*x + 5693.2*x^2
        + 146.98*x^3 + 128.19*x^4 − L0*L1*(1300.6 + 71.23*L1) + 543.8*x*L0^3 + 4.4136*L0
        − 0.71252*L0^2 − 126.38*L0^3 − 30.061*L0^4 + 5803.7*x^(-1) + 4776.5*x^(-1)*L0
        + 1001.89*x^(-1)*L0^2 + 3712/3*x^(-1)*L0^3 + 256*x^(-1)*L0^4 
        + nf*(
          80/81*L1^3 + 1040/81*L1^2
        − 16.914*L1 − 871.3 + 790.13*x − 241.23*x^2 + 43.252*x^3 − 4.3465*x*L0^3
        + 55.048*L0*L1 − 492*L0 − 343.1*L0^2 − 48.60*L0^3 + 6.0041*x^(-1) + 141.93*x^(-1)*L0
        + 2912/27*x^(-1)*L0^2 + 1280/81*x^(-1)*L0^3
        )
    )

    return value
end

function Pgq_integrand(; x::Float64, α::Float64, order::Int64, nf::Int64, f::Function)

    Pgq = Pgq0_func(x,nf)
    integrand0 = (
        Pgq*f(x)
    )

    Pgq = Pgq1_func(x,nf)
    integrand1 = (
        Pgq*f(x)
    )

    Pgq = Pgq2_func(x,nf)
    integrand2 = (
        Pgq*f(x)
    )

    #if order == 0
    #    total = (α/(4*π))*integrand0
    #elseif order == 1
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1
    #elseif order == 2
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1 + (α/(4*π))^3*integrand2
    #end

    if order == 0
        total = integrand0
    elseif order == 1
        total = integrand1
    elseif order == 2
        total = integrand2
    end
    return total
end

#Pgg-----------------------------------------------------------------------------------------------

function Pgg0_Reg_func(x,nf)

    value = -24 + 12/x + 12*x - 12*x^2 

    return value
end

function Pgg0_Delta_func(x,nf)

    value = 11 - (2*nf)/3

    return value
end

function Pgg0_D0_func(x,nf)

    value = 12

    return value
end

function Pgg1_Reg_func(x,nf)

    H0 = H0_func(x)
    H2 = H2_func(x)
    Hm10 = Hm10_func(x)
    H00 = H00_func(x)   
    H10 = H10_func(x)

    value = (-4*nf*(1 + x)*(23 + x*(-189 + x*(-45 + 121*x))) + 54*x*(-((1 + x)*(25 + 109*x)) + 6*π^2*(3 + 2*x*(2 + x + x^2))) + 
    12*(1 + x)*(-27*(22 + x*(33 + x*(3 + 22*x))) + 2*nf*(-2 + x*(57 + x*(15 + 34*x))))*H0 - 3888*H10 - 3888*H2 + 3888*Hm10 + 
    72*((-108 + x*(-81 + 54*(-4 + x)*x*(1 + x) + 4*nf*(1 + x)^2))*H00 + 54*x*((1 + x + x^3)*(H10 + H2) + (1 + x)*(2 + x + x^2)*Hm10)))/(27*x*(1 + x))

    return value
end

function Pgg1_Delta_func(x,nf)

    value = 96 - (32*nf)/3 + 108*z3

    return value
end

function Pgg1_D0_func(x,nf)

    H0 = H0_func(x)
    H2 = H2_func(x)
    H00 = H00_func(x)   
    H10 = H10_func(x)

    value = 8*(33 - 2*nf)*H0 - (4*(-201 + 10*nf + 9*π^2 + 162*H00 + 108*H10 + 108*H2))/3

    return value
end

function Pgg2_Reg_func(x,nf)

    L0 = log(x)
    L1 = log(1-x)

    value = (
        14214.4 - 804.13*nf + 1.94238*nf^2 - 28489*x + 248.949*nf*x - 77.1899*nf^2*x + 7469*x^2 + 260.6*nf*x^2 + 153.27*nf^2*x^2 + 30421*x^3 + 
        272.79*nf*x^3 - 106.029*nf^2*x^3 - 53017*x^4 + 2133.2*nf*x^4 + 11.995*nf^2*x^4 + 19556*x^5 - 926.869*nf*x^5 + 
        (3168 + 3281.7*x - 5.037*nf^2*x + 5685.8*x^2 + nf*(49.7777 + 155.1*x + 485.180*x^2))*L0^3 + (576 + (191.99 + 18.085*nf)*x)*L0^4 + 
        (-3590.1 + 319.97*nf)*x*L1 + L0^2*(3651.1 + nf^2*(1.18518 - 44.8*x) + 13528*x + nf*(263.111 + 482.94*x) + (-21328 - 29.7090*nf - 62.908*nf^2)*x*L1) + 
        L0*(10233. + nf^2*(4.54320 - 69.712*x) + 12258*x + nf*(-5.47 + 4.9934*x) + x*(-186.4 + 1266.5*nf + nf^2*(-115.009 + 96.5219*x))*L1 + 87.771*nf*x*L1^2)
    )/x

    return value
end

function Pgg2_Delta_func(x,nf)

    value = 4425.451 - 528.719*nf + 6.4628*nf^2

    return value
end

function Pgg2_D0_func(x,nf)

    value = 2643.52 - 412.172*nf - 1.77777*nf^2

    return value
end

function Pgg_integrand(; x::Float64, α::Float64, order::Int64, nf::Int64, f::Function)

    Pgg_Reg = Pgg0_Reg_func(x,nf)
    Pgg_Delta = Pgg0_Delta_func(x,nf)
    Pgg_D0_x = Pgg0_D0_func(x,nf)
    Pgg_D0_1 = Pgg0_D0_func(1,nf)
    integrand0 = (
          Pgg_Reg*f(x)
        + Pgg_Delta*f(1) #This works because integration interval is 1
        + 1/(1-x)*(Pgg_D0_x*f(x) - Pgg_D0_1*f(1))
    )

    Pgg_Reg = Pgg1_Reg_func(x,nf)
    Pgg_Delta = Pgg1_Delta_func(x,nf)
    Pgg_D0_x = Pgg1_D0_func(x,nf)
    Pgg_D0_1 = Pgg1_D0_func(1,nf)
    integrand1 = (
          Pgg_Reg*f(x)
        + Pgg_Delta*f(1) #This works because integration interval is 1
        + 1/(1-x)*(Pgg_D0_x*f(x) - Pgg_D0_1*f(1))
    )

    Pgg_Reg = Pgg2_Reg_func(x,nf)
    Pgg_Delta = Pgg2_Delta_func(x,nf)
    Pgg_D0_x = Pgg2_D0_func(x,nf)
    Pgg_D0_1 = Pgg2_D0_func(1,nf)
    integrand2 = (
          Pgg_Reg*f(x)
        + Pgg_Delta*f(1) #This works because integration interval is 1
        + 1/(1-x)*(Pgg_D0_x*f(x) - Pgg_D0_1*f(1))
    )

    #if order == 0
    #    total = (α/(4*π))*integrand0
    #elseif order == 1
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1
    #elseif order == 2
    #    total = (α/(4*π))*integrand0 + (α/(4*π))^2*integrand1 + (α/(4*π))^3*integrand2
    #end

    if order == 0
        total = integrand0
    elseif order == 1
        total = integrand1
    elseif order == 2
        total = integrand2
    end

    return total
end

# Getting the moments
using QuadGK

function Pqq_moment(; type::String, channel::String, order::Integer, nf::Integer)
    α = 1.0
    eps = 1e-12

    if type == "d0"
        f = x -> x^2
    elseif type == "d1"
        f = x -> x^2*log(x)
    elseif type == "d2"
        f = x -> x^2*log(x)^2
    else
        throw(ArgumentError("type must be one of \"d0\", \"d1\", or \"d2\""))
    end

    if channel == "qq"
        kernel = Pqq_integrand
    elseif channel == "qg"
        kernel = Pqg_integrand
    elseif channel == "gq"
        kernel = Pgq_integrand
    elseif channel == "gg"
        kernel = Pgg_integrand
    else
        throw(ArgumentError("channel must be one of \"qq\", \"qg\", \"gq\", or \"gg\""))
    end

    integrand(x) = kernel(x = x, α = α, order = order, nf = nf, f = f)
    return quadgk(integrand, Float64(eps), 1.0 - Float64(eps))
end
