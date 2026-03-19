using QCDNUM
using OffsetArrays
using SpecialFunctions
include("LIKEn.jl")
include("QCDNUM setting.jl")

#
#const A = 12.0
μ2=4.0

const gamma = 2.1953732304972053     
const g3f   = 0.21999864082553472     
const g3D   = 1.9005551041957672E-002
const Nq1   = 0.25573866909683662     
const gq1   = 6.3511479101485691E-003
const dq1   = 0.18404030632272672
const Nq2   = 0.15627533798111459     
const gq2   = 1.1501650731427362
const dq2   = 0.47402585933913333 
const p_10  = 0.0
const p_11  = 0.0
const p_12  = 0.0 

# Initialize QCDNUM -------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------------------------

#pdf = QCDNUM.allfxq(QCD_itype, z, μ2, 0, 1) 
#T = pdf[1]/z
#B = pdf[2]/z
#C = pdf[3]/z
#S = pdf[4]/z
#D = pdf[5]/z
#U = pdf[6]/z
#GL = pdf[7]/z
#UB = pdf[8]/z
#DB = pdf[9]/z
#SB = pdf[10]/z
#CB = pdf[11]/z
#BB = pdf[12]/z
#TB = pdf[13]/z

#print(U," ", UB," ", D," ", DB," ", S," ", SB," ", C," ", CB," ", B," ", BB," ", GL)