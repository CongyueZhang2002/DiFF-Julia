using OffsetArrays
const fortran_pdfs = joinpath(@__DIR__, "pdfs.so")

function set_lhapdf(i_set::Int64, name::String, i_member::Int64)
    ## i_set   : set number for LHAPDF
    ## name    : LHAPDF set name
    ## i_member: LHAPDF set member
    _name = rpad(name, 120) ## pad with spaces if needed
    ccall((:set_lhapdf_, fortran_pdfs), Cvoid,
          (Ref{Int32}, Cstring, Ref{Int32}),
          Int32(i_set), _name, Int32(i_member))
end

function get_pdfs(i_set::Int64, x::Float64, Q::Float64)
    results_base = zeros(Float64, 11)
    ccall((:get_pdfs_, fortran_pdfs), Cvoid,
          (Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}),
          Int32(i_set), x, Q, results_base)
    results = OffsetArray(results_base, -5:5)

    return results
end

function get_alphas(i_set::Int64, Q::Float64)
    alphas = Ref{Float64}()
    ccall((:get_alphas_, fortran_pdfs), Cvoid,
          (Ref{Int32}, Ref{Float64}, Ref{Float64}),
          Int32(i_set), Q, alphas)

    return alphas[]
end

function get_heavy_masses(i_set::Integer)
    mc = Ref{Cdouble}()
    mb = Ref{Cdouble}()
    mt = Ref{Cdouble}()

    ccall((:get_heavy_masses_, fortran_pdfs), Cvoid,
          (Ref{Int32}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
          Int32(i_set), mc, mb, mt)

    return (mc[], mb[], mt[])
end

