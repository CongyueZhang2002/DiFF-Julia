#----------------------------------------------------------------------------
# PDF Set
#----------------------------------------------------------------------------

const pdf_dict_0 = Dict(
    "pdfset_name" => "CT18NLO",
    "iset" => 0,
    "i_member" => 0 
) # also for strong coupling constant

const pdf_dict_1 = Dict(
    "pdfset_name" => "JAMDiFF23-transversity_lo",
    "iset" => 1,
    "i_member" => 0 
)

wdir = "results/wLQCD"
results_dir_name = "AUT_EEC_EIC"
scan_grid_name = "a_b_uniform_area_100"
dict_raw_DiFF_name = "dict_raw_DiFF"
JAM_DiFF_extrapolation_policy = :zero # :zero, :warn_zero, or :error

const pdf_dict_array = [pdf_dict_0, pdf_dict_1]

#----------------------------------------------------------------------------
# masses
#----------------------------------------------------------------------------

# quarks
const mc = 1.27
const mb = 4.18
const mt = 172.0

# bosons
const mz = 91.1876

# hadrons
const mpi = 0.13957
