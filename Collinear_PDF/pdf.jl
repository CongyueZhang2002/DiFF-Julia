include("lhapdf.jl")

using OffsetArrays

for pdf_dict in pdf_dict_array
    set_lhapdf(pdf_dict["iset"], pdf_dict["pdfset_name"], pdf_dict["i_member"])
end

function PDF_func(; iset::Int, x::Float64, μ::Float64, flavor::String)

    #if μ < 1.3
    #    μ_safe = 1.3
    #else
        μ_safe = μ 
    #end

    if flavor == "tb"
        return 0.0
    elseif flavor == "bb"
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
    elseif flavor == "t"
        return 0.0
    else
        error("Invalid flavor: $flavor")
    end

    value = get_pdfs(iset, x, μ_safe)[i]

    return value
end