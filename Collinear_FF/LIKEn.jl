#vacuum FF: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.91.014035

using OffsetArrays
using SpecialFunctions

#const A = 12.0

#const gamma = 2.1953732304972053     
#const g3f   = 0.21999864082553472     
#const g3D   = 1.9005551041957672E-002
#const Nq1   = 0.25573866909683662     
#const gq1   = 6.3511479101485691E-003
#const dq1   = 0.18404030632272672
#const Nq2   = 0.15627533798111459     
#const gq2   = 1.1501650731427362
#const dq2   = 0.47402585933913333 
#const p_10  = 0.0
#const p_11  = 0.0
#const p_12  = 0.0 

function nFF_func(ipdf,z)::Float64

    Ni = OffsetArray(zeros(7), 0:6)
    ai = OffsetArray(zeros(7), 0:6)
    bi = OffsetArray(zeros(7), 0:6)
    gi = OffsetArray(zeros(7), 0:6)
    di = OffsetArray(zeros(7), 0:6)

    #DEHSS-----------------------------------------------------------------------------------------

    #utot
    Ni[1]  = 0.387 
    ai[1]  = -0.388
    bi[1]  = 0.910
    gi[1]  = 7.15
    di[1]  = 3.96

    #dtot
    Ni[2]  = 0.388 
    ai[2]  = ai[1]
    bi[2]  = bi[1]
    gi[2]  = gi[1]
    di[2]  = di[1]

    #ubar = d
    Ni[3] = 0.105 
    ai[3] = 1.649
    bi[3] = 3.286
    gi[3] = 49.95
    di[3] =  8.67

    #stot
    Ni[4]  = 0.273 
    ai[4]  = 1.449
    bi[4]  = bi[3]
    gi[4]  = gi[3]
    di[4]  = di[3]

    #ctot
    Ni[5]  = 0.306 
    ai[5]  = 1.345
    bi[5]  = 5.519
    gi[5]  = 19.78
    di[5]  = 10.22

    #btot
    Ni[6]  =  0.372 
    ai[6]  = -0.127 
    bi[6]  =  4.490 
    gi[6]  =  24.49 
    di[6]  =  12.80

    #g
    Ni[0]   =  0.260 
    ai[0]   =  2.552
    bi[0]   =  6.194
    gi[0]   =  87.06
    di[0]   =  20.36  

    #Nuclear Modification--------------------------------------------------------------------------

    #quarks
    for i in 1:6
        Ni[i] = Ni[i] / (beta(2.0+ai[i], bi[i]+1.0) + gi[i] * beta(2.0+ai[i], bi[i]+di[i]+1.0)) #maybe should be put later   
        Ni[i] = Ni[i] * (1.0 + Nq1*(1.0-A^Nq2))
        ai[i] = ai[i] 
        bi[i] = bi[i] 
        gi[i] = gi[i] + gq1 *(1.0-A^gq2)
        di[i] = di[i] + dq1 *(1.0-A^dq2)     
    end


    #----------------------------------------------------------------------------------------------

    #gluon
    Ni[0] = Ni[0] / (beta(2.0+ai[0], bi[0]+1.0) + gi[0] * beta(2.0+ai[0], bi[0]+di[0]+1.0))
    Ni[0] = Ni[0] 
    ai[0] = ai[0] 
    bi[0] = bi[0] 
    gi[0] = gi[0] 
    di[0] = di[0] 
          
    utot  = Ni[1]*z^ai[1]*(1-z)^bi[1]*(1.0+gi[1]*(1-z)^di[1])
    dtot  = Ni[2]*z^ai[2]*(1-z)^bi[2]*(1.0+gi[2]*(1-z)^di[2])
    ub    = Ni[3]*z^ai[3]*(1-z)^bi[3]*(1.0+gi[3]*(1-z)^di[3]) 
    d     = ub                                                
    stot  = Ni[4]*z^ai[4]*(1-z)^bi[4]*(1.0+gi[4]*(1-z)^di[4])
    ctot  = Ni[5]*z^ai[5]*(1-z)^bi[5]*(1.0+gi[5]*(1-z)^di[5])
    btot  = Ni[6]*z^ai[6]*(1-z)^bi[6]*(1.0+gi[6]*(1-z)^di[6])
    g     = Ni[0]*z^ai[0]*(1-z)^bi[0]*(1.0+gi[0]*(1-z)^di[0]) 

    u   = utot - ub 
    db  = dtot - d  
    s   = stot/2.0  
    sb  = s         
    c   = ctot/2.0  
    cb  = c         
    b   = btot/2.0
    bb  = b

    func=0.0
    if(ipdf == 0) 
        func = g
    elseif(ipdf == 1) 
        func = d
    elseif(ipdf == 2) 
        func = u
    elseif(ipdf == 3) 
        func = s 
    elseif(ipdf == 4) 
        func = db 
    elseif(ipdf == 5) 
        func = ub 
    elseif(ipdf == 6) 
        func = sb   
    elseif(ipdf == 7) 
        func = c
    elseif(ipdf == 8) 
        func = cb 
    elseif(ipdf == 9) 
        func = b 
    elseif(ipdf == 10) 
        func = bb  
    elseif(ipdf == 11) 
        func = 0.0
    elseif(ipdf == 12) 
        func = 0.0
    end
    return z*func
end

#print(nFF_func(1,0.1))

function FF_deuteron_func(; z::Float64, μ::Float64, flavor::String)

    if μ < sqrt(1.005)
        μ_safe = 1.005
    else
        μ_safe = μ 
    end

    if flavor == "bb"
        i = -5
    elseif flavor == "cb"
        i = -4
    elseif flavor == "sb"
        i = -3
    elseif flavor == "db"
        i = -2
    elseif flavor == "ub"
        i = -1
    elseif flavor == "g"
        i = 0
    elseif flavor == "u"
        i = 1
    elseif flavor == "d"
        i = 2
    elseif flavor == "s"
        i = 3
    elseif flavor == "c"
        i = 4
    elseif flavor == "b"
        i = 5
    end

    value = get_pdfs(-2, z, μ_safe)[i]

    return value
end