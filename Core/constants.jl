# masses

#const mc, mb, mt = get_heavy_masses(iset)

# TMD

const b0 = 1.12292

# SU(3)

const NC = 3
const CF = 4/3
const CA = 3
const TF = 1/2

const NA = NC^2 - 1 
const NF = NC 

# SU(3) Casimir Tensors

const dAA_NA = 135/8
const dFA_NA = 15/16
const dFF_NA = 5/96

const dFA_NF = dFA_NA * NA / NF 
const dFF_NF = dFF_NA * NA / NF 

# SU(3) four-loop colour-coefficients from arXiv:2205.02249 (Table 1)

const b_d4_FA = 998.0       #  b(d^{(4)}_{FA})
const b_d4_FF = 143.6       #  b(d^{(4)}_{FF})
const b_CA_CF2_nf = 455.247 #  b(N_f * C_A * C_F^2)
const b_CF3_nf = -80.780    #  b(N_f * C_F^3)

# Strong coupling
const ΛQCD3 = 0.326
const ΛQCD4 = 0.326
const ΛQCD5 = 0.226

# Quarks
const eu = 2/3
const ed = 1/3

#EM
const αem = 1/137

# Weak
const sw = sqrt(0.2313)
const cw = sqrt(1-sw*sw)

const GZ = 2.5

# Riemann Zeta Constants

const z2 = 1.644934067
const z3 = 1.202056903
const z4 = 1.082323234
const z5 = 1.036927755
const z6 = 1.017343062
const z7 = 1.008349277

# Numerical

const rtol = 1e-3
const DY_conversion = 0.3894*10^9
const b_thres = 30.0

rtol_TMD = 5e-3
rtol_Sudakov = 1e-3
