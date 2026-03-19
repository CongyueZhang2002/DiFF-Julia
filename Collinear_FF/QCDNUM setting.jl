using QCDNUM

# Set-up 

QCD_itype = 3 # 1: Unpol PDF, 2: Pol PDF, 3: FF 
QCD_iord = 2 # 1 : LO, 2 : NLO, 3 : NNLO
QCD_iosp = 2 # 2 for linear interpolarion, 3 for spline
QCD_q0 = 1.0 # the scale of the initial nFF 

# Physics parameters

QCD_r20 = 2.0 #
QCD_as0 = 0.364 # value of strong coupling at scale r20
QCD_q2c = (1.43)^2 # charm quark mass squared
QCD_q2b = (4.3)^2 # bottom quark mass squared

# x grid

QCD_xmin = Float64.([1.0e-4]) 
QCD_nxin = 100
QCD_iwt = Int32.([1])
QCD_ng = 1

# μ grid    

QCD_qq = Float64.([1e0, 1e5]) 
QCD_wt = Float64.([1e0, 1e0])
QCD_nqin = 60 
QCD_ngq = 2

QCD_iqt = 999
QCD_nfin = 1

#              tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
QCD_def = Float64.([0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., # d  
                    0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., # u      
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., # s      
                    0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., # db     
                    0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., # ub    
                    0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., # sb     
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., # c     
                    0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # cb     
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., # b     
                    0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # bb     
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., # t     
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.] # tb
);

QCD_wrapped_func = WrappedPDF(nFF_func)

QCDNUM_initialization_sequence = """
QCDNUM.qcinit(-6, " ")
nx = QCDNUM.gxmake(QCD_xmin, QCD_iwt, QCD_ng, QCD_nxin, QCD_iosp)
nq = QCDNUM.gqmake(QCD_qq, QCD_wt, QCD_ngq, QCD_nqin)
nw = QCDNUM.fillwt(QCD_itype)
QCDNUM.setord(QCD_iord)
QCDNUM.setalf(QCD_as0, QCD_r20)
QCD_iqc = QCDNUM.iqfrmq(QCD_q2c)
QCD_iqb = QCDNUM.iqfrmq(QCD_q2b)
QCDNUM.setcbt(QCD_nfin, QCD_iqc, QCD_iqb, QCD_iqt)
QCD_iq0 = QCDNUM.iqfrmq(QCD_q0)
eps = QCDNUM.evolfg(QCD_itype, QCD_wrapped_func, QCD_def, QCD_iq0)
"""