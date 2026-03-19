subroutine set_lhapdf(set, pdf_name, member)
    character*120 pdf_name
    integer :: set, member

    ! call SetLHAPARM('silent') !! suppress messages from LHAPDF
    call InitPDFsetByNameM(set, trim(pdf_name))
    call InitPDFM(set, member)
end subroutine set_lhapdf

subroutine get_pdfs(set, x, Q, results)
    integer set
    double precision x, Q
    double precision xPDFs(-6:6), results(-5:5)

    call evolvePDFM(set, x, Q, xPDFs)

    results( 0) = xPDFs( 0) / x
    results( 1) = xPDFs( 2) / x
    results( 2) = xPDFs( 1) / x
    results( 3) = xPDFs( 3) / x
    results( 4) = xPDFs( 4) / x
    results( 5) = xPDFs( 5) / x
    results(-1) = xPDFs(-2) / x
    results(-2) = xPDFs(-1) / x
    results(-3) = xPDFs(-3) / x
    results(-4) = xPDFs(-4) / x
    results(-5) = xPDFs(-5) / x
    return
end subroutine get_pdfs

subroutine get_alphas(i_set, Q, alphas)
    !! returns alpha_s
    implicit none
    integer i_set
    double precision Q, alphas, alphasPDFM
    external alphasPDF

    alphas = alphasPDFM(i_set, Q)
    return
end subroutine get_alphas

subroutine get_heavy_masses(i_set, mc, mb, mt)
    !! returns (mc, mb, mt) from the PDF set i_set
    implicit none
    integer           :: i_set
    double precision  :: mc, mb, mt

    call GetQmassM(i_set, 4, mc)   ! charm
    call GetQmassM(i_set, 5, mb)   ! bottom
    call GetQmassM(i_set, 6, mt)   ! top

    return
end subroutine get_heavy_masses