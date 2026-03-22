#qq_x^2--------------------------------------------------------------------------------------------

function γqq0_func(nf)

    if nf == 3
        value = 5.555556
    elseif nf == 4
        value = 5.555556
    elseif nf == 5
        value = 5.555556
    elseif nf == 6
        value = 5.555556
    end

    return value
end

function γqq1_func(nf)

    if nf == 3
        value = 40.517246
    elseif nf == 4
        value = 36.477987
    elseif nf == 5
        value = 32.438727
    elseif nf == 6
        value = 28.399468
    end

    return value
end

function γqq2_func(nf)

    if nf == 3
        value = 372.35148
    elseif nf == 4
        value = 212.183306
    elseif nf == 5
        value = 50.273165
    elseif nf == 6
        value = -113.378943
    end

    return value
end

#qg_x^2--------------------------------------------------------------------------------------------

function γqg0_func(nf)

    if nf == 3
        value = -1.4
    elseif nf == 4
        value = -1.866667
    elseif nf == 5
        value = -2.333333
    elseif nf == 6
        value = -2.8
    end

    return value
end

function γqg1_func(nf)

    if nf == 3
        value = 12.8323
    elseif nf == 4
        value = 16.469733
    elseif nf == 5
        value = 19.787166
    elseif nf == 6
        value = 22.7846
    end

    return value
end

function γqg2_func(nf)

    if nf == 3
        value = -94.626677
    elseif nf == 4
        value = -114.708478
    elseif nf == 5
        value = -130.588951
    elseif nf == 6
        value = -143.185426
    end

    return value
end

#gq_x^2--------------------------------------------------------------------------------------------

function γgq0_func(nf)

    if nf == 3
        value = -1.555556
    elseif nf == 4
        value = -1.555556
    elseif nf == 5
        value = -1.555556
    elseif nf == 6
        value = -1.555556
    end

    return value
end

function γgq1_func(nf)

    if nf == 3
        value = -44.885707
    elseif nf == 4
        value = -44.885707
    elseif nf == 5
        value = -44.885707
    elseif nf == 6
        value = -44.885707
    end

    return value
end

function γgq2_func(nf)

    if nf == 3
        value = -535.840112
    elseif nf == 4
        value = -454.706098
    elseif nf == 5
        value = -373.572083
    elseif nf == 6
        value = -292.438068
    end

    return value
end

#gg_x^2--------------------------------------------------------------------------------------------

function γgg0_func(nf)

    if nf == 3
        value = 10.4
    elseif nf == 4
        value = 11.066667
    elseif nf == 5
        value = 11.733333
    elseif nf == 6
        value = 12.4
    end

    return value
end

function γgg1_func(nf)

    if nf == 3
        value = -50.545952
    elseif nf == 4
        value = -65.243663
    elseif nf == 5
        value = -79.941373
    elseif nf == 6
        value = -94.639084
    end

    return value
end

function γgg2_func(nf)

    if nf == 3
        value = -382.155909
    elseif nf == 4
        value = -523.942484
    elseif nf == 5
        value = -675.768406
    elseif nf == 6
        value = -837.633674
    end

    return value
end

function γ_func(; type::String, α::Float64, order::Int64, nf::Int64)

    if     type == "qq"
        order0 = γqq0_func(nf)
        order1 = γqq1_func(nf)
        order2 = γqq2_func(nf)
    elseif type == "qg"
        order0 = γqg0_func(nf)
        order1 = γqg1_func(nf)
        order2 = γqg2_func(nf)
    elseif type == "gq"
        order0 = γgq0_func(nf)
        order1 = γgq1_func(nf)
        order2 = γgq2_func(nf)
    elseif type == "gg"
        order0 = γgg0_func(nf)
        order1 = γgg1_func(nf)
        order2 = γgg2_func(nf)
    else
        error("Type $type not supported.")
    end

    as = α/(4*π)

    if order < 0
        total = 0
    elseif order == 0
        total = as*order0
    elseif order == 1
        total = as*order0 + as^2*order1
    elseif order == 2
        total = as*order0 + as^2*order1 + as^3*order2
    end

    return total
end

#qq_x^2*log(x)-------------------------------------------------------------------------------------

function dγqq0_func(nf)

    if nf == 3
        value = 1.643352
    elseif nf == 4
        value = 1.643352
    elseif nf == 5
        value = 1.643352
    elseif nf == 6
        value = 1.643352
    end

    return value
end

function dγqq1_func(nf)

    if nf == 3
        value = 11.222329
    elseif nf == 4
        value = 8.764223
    elseif nf == 5
        value = 6.306117
    elseif nf == 6
        value = 3.848011
    end

    return value
end

function dγqq2_func(nf)

    if nf == 3
        value = 149.52916
    elseif nf == 4
        value = 88.321096
    elseif nf == 5
        value = 25.684808
    elseif nf == 6
        value = -38.379703
    end

    return value
end

#qg_x^2*log(x)-------------------------------------------------------------------------------------

function dγqg0_func(nf)

    if nf == 3
        value = 0.396667
    elseif nf == 4
        value = 0.528889
    elseif nf == 5
        value = 0.661111
    elseif nf == 6
        value = 0.793333
    end

    return value
end

function dγqg1_func(nf)

    if nf == 3
        value = -5.889708
    elseif nf == 4
        value = -8.30236
    elseif nf == 5
        value = -10.939721
    elseif nf == 6
        value = -13.801789
    end

    return value
end

function dγqg2_func(nf)

    if nf == 3
        value = -13.565441
    elseif nf == 4
        value = -18.929944
    elseif nf == 5
        value = -25.759532
    elseif nf == 6
        value = -34.680448
    end

    return value
end

#gq_x^2*log(x)-------------------------------------------------------------------------------------

function dγgq0_func(nf)

    if nf == 3
        value = 0.907407
    elseif nf == 4
        value = 0.907407
    elseif nf == 5
        value = 0.907407
    elseif nf == 6
        value = 0.907407
    end

    return value
end

function dγgq1_func(nf)

    if nf == 3
        value = 15.135048
    elseif nf == 4
        value = 15.135048
    elseif nf == 5
        value = 15.135048
    elseif nf == 6
        value = 15.135048
    end

    return value
end

function dγgq2_func(nf)

    if nf == 3
        value = -199.25068
    elseif nf == 4
        value = -238.53264
    elseif nf == 5
        value = -277.814599
    elseif nf == 6
        value = -317.096559
    end

    return value
end

#gg_x^2*log(x)-------------------------------------------------------------------------------------

function dγgg0_func(nf)

    if nf == 3
        value = 5.342542
    elseif nf == 4
        value = 5.342542
    elseif nf == 5
        value = 5.342542
    elseif nf == 6
        value = 5.342542
    end

    return value
end

function dγgg1_func(nf)

    if nf == 3
        value = 40.809479
    elseif nf == 4
        value = 36.279598
    elseif nf == 5
        value = 31.749718
    elseif nf == 6
        value = 27.219838
    end

    return value
end

function dγgg2_func(nf)

    if nf == 3
        value = -270.251875
    elseif nf == 4
        value = -460.409564
    elseif nf == 5
        value = -655.315333
    elseif nf == 6
        value = -854.969181
    end

    return value
end

function dγ_func(; type::String, α::Float64, order::Int64, nf::Int64)

    if     type == "qq"
        order0 = dγqq0_func(nf)
        order1 = dγqq1_func(nf)
        order2 = dγqq2_func(nf)
    elseif type == "qg"
        order0 = dγqg0_func(nf)
        order1 = dγqg1_func(nf)
        order2 = dγqg2_func(nf)
    elseif type == "gq"
        order0 = dγgq0_func(nf)
        order1 = dγgq1_func(nf)
        order2 = dγgq2_func(nf)
    elseif type == "gg"
        order0 = dγgg0_func(nf)
        order1 = dγgg1_func(nf)
        order2 = dγgg2_func(nf)
    else
        error("Type $type not supported.")
    end

    as = α/(4*π)

    if order < 0
        total = 0
    elseif order == 0
        total = as*order0
    elseif order == 1
        total = as*order0 + as^2*order1
    elseif order == 2
        total = as*order0 + as^2*order1 + as^3*order2
    end

    return total
end

#qq_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγqq0_func(nf)

    if nf == 3
        value = -0.541076
    elseif nf == 4
        value = -0.541076
    elseif nf == 5
        value = -0.541076
    elseif nf == 6
        value = -0.541076
    end

    return value
end

function ddγqq1_func(nf)

    if nf == 3
        value = -0.75182
    elseif nf == 4
        value = 0.789334
    elseif nf == 5
        value = 2.330487
    elseif nf == 6
        value = 3.871641
    end

    return value
end

function ddγqq2_func(nf)

    if nf == 3
        value = -57.94726
    elseif nf == 4
        value = -46.872315
    elseif nf == 5
        value = -35.206391
    elseif nf == 6
        value = -22.949487
    end

    return value
end

#qg_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγqg0_func(nf)

    if nf == 3
        value = -0.261444
    elseif nf == 4
        value = -0.348593
    elseif nf == 5
        value = -0.435741
    elseif nf == 6
        value = -0.522889
    end

    return value
end

function ddγqg1_func(nf)

    if nf == 3
        value = 5.948939
    elseif nf == 4
        value = 8.254236
    elseif nf == 5
        value = 10.720691
    elseif nf == 6
        value = 13.348304
    end

    return value
end

function ddγqg2_func(nf)

    if nf == 3
        value = -13.438896
    elseif nf == 4
        value = -25.170153
    elseif nf == 5
        value = -39.80025
    elseif nf == 6
        value = -56.893
    end

    return value
end

#gq_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγgq0_func(nf)

    if nf == 3
        value = -1.021605
    elseif nf == 4
        value = -1.021605
    elseif nf == 5
        value = -1.021605
    elseif nf == 6
        value = -1.021605
    end

    return value
end

function ddγgq1_func(nf)

    if nf == 3
        value = -0.135984
    elseif nf == 4
        value = -0.135984
    elseif nf == 5
        value = -0.135984
    elseif nf == 6
        value = -0.135984
    end

    return value
end

function ddγgq2_func(nf)

    if nf == 3
        value = 717.599853
    elseif nf == 4
        value = 747.709978
    elseif nf == 5
        value = 777.820103
    elseif nf == 6
        value = 807.930228
    end

    return value
end

#gg_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγgg0_func(nf)

    if nf == 3
        value = -3.254588
    elseif nf == 4
        value = -3.254588
    elseif nf == 5
        value = -3.254588
    elseif nf == 6
        value = -3.254588
    end

    return value
end

function ddγgg1_func(nf)

    if nf == 3
        value = 5.256067
    elseif nf == 4
        value = 7.301451
    elseif nf == 5
        value = 9.346834
    elseif nf == 6
        value = 11.392218
    end

    return value
end

function ddγgg2_func(nf)

    if nf == 3
        value = 1288.956232
    elseif nf == 4
        value = 1378.800628
    elseif nf == 5
        value = 1472.630253
    elseif nf == 6
        value = 1570.445109
    end

    return value
end

function ddγ_func(; type::String, α::Float64, order::Int64, nf::Int64)

    if     type == "qq"
        order0 = ddγqq0_func(nf)
        order1 = ddγqq1_func(nf)
        order2 = ddγqq2_func(nf)
    elseif type == "qg"
        order0 = ddγqg0_func(nf)
        order1 = ddγqg1_func(nf)
        order2 = ddγqg2_func(nf)
    elseif type == "gq"
        order0 = ddγgq0_func(nf)
        order1 = ddγgq1_func(nf)
        order2 = ddγgq2_func(nf)
    elseif type == "gg"
        order0 = ddγgg0_func(nf)
        order1 = ddγgg1_func(nf)
        order2 = ddγgg2_func(nf)
    else
        error("Type $type not supported.")
    end

    as = α/(4*π)

    if order < 0
        total = 0
    elseif order == 0
        total = as*order0
    elseif order == 1
        total = as*order0 + as^2*order1
    elseif order == 2
        total = as*order0 + as^2*order1 + as^3*order2
    end

    return total
end
