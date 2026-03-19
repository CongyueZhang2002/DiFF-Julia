using Distributed
using Printf
using Statistics

#----------------------------------------------------------------------------
# Initialize Workers
#----------------------------------------------------------------------------

N_threads = 16
addprocs(N_threads)

@everywhere begin
#----------------------------------------------------------------------------
# Which Fit?
#----------------------------------------------------------------------------

fit_name = "MSHT20N3LO-MC"
if_grid = false

#----------------------------------------------------------------------------
# Include Files
#----------------------------------------------------------------------------

inc(path) = include(joinpath(@__DIR__, path))
# Fit Card
inc("../Cards/$fit_name.jl")
# grid and interpolations
inc("grid.jl")
inc("interpolation.jl")
# PDF
inc("../Collinear_PDF/pdf.jl")
# Core
inc("../Core/constants.jl")
inc("../Core/strong coupling.jl")
# Numerical
inc("../Numerical/FastGK.jl")
# TMDs
inc("../TMDs/NP Parameterizations/$NP_name")
inc("../TMDs/TMDPDFs/$TMDPDF_name")

end

#----------------------------------------------------------------------------
# Settings for grid generation
#----------------------------------------------------------------------------

# independent variables

variable_names = ["x", "b"]

eps_x = 1e-5
eps_b = 1e-4
variable_ranges = [
    [eps_x,1-eps_x],
    [eps_b,25.0],
]
 
variable_settings = [
    [300, "power", exp(1), 3.0],
    [400, "power", exp(1), 4.0],
]

# outputs 

generate = true
test = false

target_names = ["f_u", "f_ub", "f_d", "f_db", "f_s", "f_sb" 
                , "f_c", "f_cb", "f_b", "f_bb"]

@everywhere target_func(x, b) = xTMDPDF_raw_func(x, b)

precision = 6

path_save = "$fit_name/TMDPDF.csv"
path = joinpath(@__DIR__, path_save)
mkpath(dirname(path))  
isfile(path) && rm(path)

# testing grid

interpolator_name = "example"

N_samples = 1000

threshold = 1e-4

#----------------------------------------------------------------------------
# Run
#----------------------------------------------------------------------------

if generate == true
    grid_generator(
        variable_ranges=variable_ranges, 
        variable_settings=variable_settings, 
        variable_names=variable_names,                
        target_func=target_func, 
        target_names=target_names, 
        precision=precision, 
        path=path
    )
end

#----------------------------------------------------------------------------
# Testing
#----------------------------------------------------------------------------

if test == true

    function random_point_generator(; variable_ranges, variable_settings, N_samples)

        N_dim = length(variable_ranges)

        arrays = [
            array_generator(
                a0    = float(variable_ranges[d][1]),
                an    = float(variable_ranges[d][2]),
                n     = variable_settings[d][1],
                method= variable_settings[d][2],
                base  = variable_settings[d][3],
                power = variable_settings[d][4]
            ) for d in 1:N_dim
        ]

        rand_in_bin(arr) = begin
            idx = rand(1:length(arr)-1)   
            lo, hi = arr[idx], arr[idx+1]
            lo + (hi - lo) * rand()      
        end

        Point = NTuple{N_dim, Float64}
        out   = Vector{Point}(undef, N_samples)
        for i in 1:N_samples
            out[i] = ntuple(d -> rand_in_bin(arrays[d]), N_dim)
        end

        return out
    end

    df = DataFrame(CSV.File(path))
    Main.initialize_interpolator(df=df, interpolator_name=interpolator_name, 
                                variable_names=variable_names, target_names=target_names)

    points = random_point_generator(variable_ranges=variable_ranges, variable_settings=variable_settings, N_samples=N_samples)

    points_tuples = points isa AbstractMatrix ? [Tuple(row) for row in eachrow(points)] : collect(points)
    N = length(points_tuples)

    exact_values        = Dict(name => Vector{Float64}(undef, N) for name in target_names)
    interpolated_values = Dict(name => Vector{Float64}(undef, N) for name in target_names)
    rel_errors          = Dict(name => Vector{Float64}(undef, N) for name in target_names)

    interp = interpolators[Symbol(interpolator_name)]

    print("\x1b[2J\x1b[H"); flush(stdout)
    p = Progress(N; dt=0.1, desc="Testing...", output=stdout)
    for (i, pt) in enumerate(points_tuples)
        y_hat  = interp(pt...)           
        y_true = target_func(pt...)      

        yh = y_hat  isa Number ? (y_hat,)  : Tuple(y_hat)
        yt = y_true isa Number ? (y_true,) : Tuple(y_true)

        for (j, name) in enumerate(target_names)
            exact_values[name][i]        = yt[j]
            interpolated_values[name][i] = yh[j]
            rel_errors[name][i]          = (yh[j] - yt[j]) / yt[j] * 100.0
        end
        next!(p)
    end
    finish!(p); println()

    points_filtered = Dict{String, Vector{typeof(points_tuples[1])}}()
    for name in target_names
        m   = abs.(exact_values[name]) .>= threshold
        idx = findall(m)
        exact_values[name]        = exact_values[name][idx]
        interpolated_values[name] = interpolated_values[name][idx]
        rel_errors[name]          = rel_errors[name][idx]
        points_filtered[name]     = [points_tuples[k] for k in idx]
    end

    r2(x) = parse(Float64, @sprintf("%.2g", x))
    r5(x) = parse(Float64, @sprintf("%.5g", x))
    qs = (99.9, 99.0, 90.0, 50.0)

    for name in target_names
        a   = abs.(Float64.(rel_errors[name]))
        pts = points_filtered[name]

        println(name, ":")
        println(length(pts), " points left")

        finite = isfinite.(a)
        if any(finite)
            inds = findall(finite); gi = inds[argmax(a[finite])]
            coords = join(string.(r5.(collect(pts[gi]))), ", ")
            y_true = exact_values[name][gi]; y_hat = interpolated_values[name][gi]
            @printf("maximum error: %.2g%% at point (%s)\n", a[gi], coords)
            @printf("  y_true=%.5g, y_hat=%.5g\n", y_true, y_hat)
        else
            println("maximum error: (no finite values)")
        end

        af = a[finite]
        for q in qs
            v = isempty(af) ? NaN : quantile(af, q/100)
            @printf("%g%% percentile error: %.2g%%\n", q, v)
        end
        println()
    end

end

rmprocs(workers())
println("workers removed")