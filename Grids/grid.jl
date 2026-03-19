using DataFrames, CSV, ProgressMeter, ProgressBars

function array_generator(; a0::Float64, an::Float64, n::Integer, method::String, 
    base::Float64=exp(1.0), power::Float64=2.0)

    if a0 > an || n < 2 || base <= 1.0 || power <= 0.0 
        error("Input Invalid")
    elseif method == "linear"
        array = range(a0, an; length=n)
    elseif method == "log"
        if a0 <= 0.0 || an <= 0.0
            error("For log method, a0 and an must be positive")
        end
        a0_lin = log(a0)/log(base) 
        an_lin = log(an)/log(base)
        array = base .^ range(a0_lin, an_lin; length=n)
    elseif method == "power"
        u = range(0.0, 1.0; length=n)
        array = a0 .+ (an - a0) .* (u .^ power)
    else
        error("Choose method from linear, log, or power; log has base option, power has power option")
    end
    return collect(array)
end

function array_assembler(; variable_ranges, variable_settings)

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

    return vec(collect(Iterators.product(arrays...)))
end

function grid_generator_raw(; variables_grid, variable_names, target_func, target_names, precision)

    N_variables  = length(first(variables_grid))
    df = DataFrame((Symbol(variable_names[i]) => getindex.(variables_grid, i) for i in 1:N_variables)...)

    print("\x1b[2J\x1b[H"); flush(stdout)
    print("Generating Grid with $(length(variables_grid)) points...\n")
    predictions = pmap(t -> target_func(t...), ProgressBar(variables_grid); batch_size=100)
    print("Grid Generation Complete!\n")

    predictions = map(v -> v isa Number ? (v,) : v, predictions)  

    for (j, nm) in enumerate(target_names)
        df[!, Symbol(nm)] = getindex.(predictions, j)
    end

    transform!(df, names(df, Number) .=> ByRow(x -> round(x; sigdigits=precision)) .=> names(df, Number))

    return df
end

function grid_generator(; variable_ranges, variable_settings, variable_names, 
                        target_func, target_names, precision=6, path)

    variables_grid = array_assembler(variable_ranges=variable_ranges, variable_settings=variable_settings)

    df = grid_generator_raw(
        variables_grid = variables_grid, 
        variable_names = variable_names, 
        target_func = target_func, 
        target_names = target_names,
        precision = precision
    )

    CSV.write(path, df)
end