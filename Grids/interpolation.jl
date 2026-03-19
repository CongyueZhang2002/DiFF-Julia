using CSV, DataFrames, Interpolations, StaticArrays

const interpolators = Dict{Symbol,Any}()

function _regular_axis(axis::AbstractVector{<:AbstractFloat})
    if length(axis) <= 1
        return range(first(axis), length = length(axis))
    end

    step = axis[2] - axis[1]
    for i in 3:length(axis)
        if !isapprox(axis[i] - axis[i - 1], step; atol = 1f-5, rtol = 1f-5)
            throw(ArgumentError("cubic interpolation requires regularly spaced axes"))
        end
    end

    return range(first(axis), step = step, length = length(axis))
end

function initialize_interpolator(; df, interpolator_name, variable_names, target_names, method::Symbol = :linear)
    # axes & values in Float32 → less memory traffic at query time
    variable_axes = map(v -> Float32.(sort(unique(df[!, v]))), variable_names)
    shape = map(length, variable_axes)
    N_targets = length(target_names)

    target_vec = Array{SVector{N_targets,Float32}}(undef, shape...)
    for row in eachrow(df)
        grid_index = ntuple(d -> searchsortedfirst(variable_axes[d], Float32(row[variable_names[d]])),
                            length(variable_names))
        target_vec[grid_index...] = SVector{N_targets,Float32}((Float32(row[t]) for t in target_names)...)
    end

    base_interpolant =
        if method == :linear
            interpolate(Tuple(variable_axes), target_vec, Gridded(Linear()))
        elseif method == :cubic
            scaled_axes = map(_regular_axis, variable_axes)
            scale(interpolate(target_vec, BSpline(Cubic(Line(OnGrid())))), scaled_axes...)
        else
            throw(ArgumentError("method must be :linear or :cubic"))
        end

    fill_value = zero(SVector{N_targets,Float32})
    vector_interpolant = extrapolate(base_interpolant, fill_value)

    interpolators[Symbol(interpolator_name)] = vector_interpolant
    return nothing
end


# Example
# evaluate_interpolator(name::Symbol, coords...) = interpolators[name](coords...)
