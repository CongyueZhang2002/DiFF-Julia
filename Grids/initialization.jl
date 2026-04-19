include("interpolation.jl")

@inline mu_key(mu::Real) = Float16(mu)

if @isdefined(dict_raw_DiFF) && @isdefined(selected_keys)
    global JAM_DiFF_interpolators = Dict{Any,Any}()
    global JAM_DiFF_bounds_warnings_seen = Set{Tuple{String,Symbol,Symbol}}()

    zvals = Float64.(dict_raw_DiFF["z_array"])
    Mhvals = Float64.(dict_raw_DiFF["Mh_array"])
    z_bounds = (minimum(zvals), maximum(zvals))
    Mh_bounds = (minimum(Mhvals), maximum(Mhvals))
    flavors_by_kind = haskey(dict_raw_DiFF, "flavors_by_kind") ? dict_raw_DiFF["flavors_by_kind"] : nothing

    function JAM_DiFF_extrapolation_policy_value()
        policy = isdefined(@__MODULE__, :JAM_DiFF_extrapolation_policy) ?
            getfield(@__MODULE__, :JAM_DiFF_extrapolation_policy) :
            :warn_zero
        policy = policy isa Symbol ? policy : Symbol(policy)
        if !(policy in (:zero, :warn_zero, :error))
            throw(ArgumentError("JAM_DiFF_extrapolation_policy must be :zero, :warn_zero, or :error"))
        end
        return policy
    end

    function check_JAM_DiFF_bounds(; kind::AbstractString, z::Real, Mh::Real, mu::Real, rep::Int)
        z_float = Float64(z)
        Mh_float = Float64(Mh)
        out_of_bounds = Tuple{Symbol,Symbol,Float64,Tuple{Float64,Float64}}[]

        if z_float < z_bounds[1]
            push!(out_of_bounds, (:z, :low, z_float, z_bounds))
        elseif z_float > z_bounds[2]
            push!(out_of_bounds, (:z, :high, z_float, z_bounds))
        end

        if Mh_float < Mh_bounds[1]
            push!(out_of_bounds, (:Mh, :low, Mh_float, Mh_bounds))
        elseif Mh_float > Mh_bounds[2]
            push!(out_of_bounds, (:Mh, :high, Mh_float, Mh_bounds))
        end

        isempty(out_of_bounds) && return nothing

        policy = JAM_DiFF_extrapolation_policy_value()
        policy == :zero && return nothing

        messages = [
            string(name, "=", value, " is ", side, "er than grid range [", bounds[1], ", ", bounds[2], "]")
            for (name, side, value, bounds) in out_of_bounds
        ]
        message = string(
            "JAM DiFF grid lookup outside interpolation support for ",
            "kind=", kind, ", mu=", mu, ", rep=", rep, ": ",
            join(messages, "; "),
            ". Set JAM_DiFF_extrapolation_policy = :zero to keep the legacy silent-zero behavior, ",
            "or regenerate the grid with wider support."
        )

        if policy == :error
            throw(ArgumentError(message))
        end

        for (name, side, _, _) in out_of_bounds
            warning_key = (String(kind), name, side)
            if !(warning_key in JAM_DiFF_bounds_warnings_seen)
                @warn message
                push!(JAM_DiFF_bounds_warnings_seen, warning_key)
            end
        end

        return nothing
    end

    function build_JAM_DiFF_interpolator!(key)
        kind = key[1]
        flavor_names = flavors_by_kind === nothing ? String.(dict_raw_DiFF["flavors"]) : String.(flavors_by_kind[kind])
        values = Float64.(dict_raw_DiFF[key])

        df = DataFrame()
        df[!, "z"] = zvals
        df[!, "Mh"] = Mhvals
        for (j, flav) in enumerate(flavor_names)
            df[!, flav] = values[:, j]
        end

        interpolator_name = "JAM_DiFF_" * replace(string(key), r"[^A-Za-z0-9]+" => "_")
        initialize_interpolator(
            df = df,
            interpolator_name = interpolator_name,
            variable_names = ["z", "Mh"],
            target_names = flavor_names,
            method = :cubic,
        )
        JAM_DiFF_interpolators[key] = interpolators[Symbol(interpolator_name)]
    end

    for key in selected_keys
        build_JAM_DiFF_interpolator!(key)
    end

    global JAM_DiFF_grid
    @inline function JAM_DiFF_grid(; kind::AbstractString, z::Real, Mh::Real, mu::Real, rep::Int)
        key = (String(kind), mu_key(mu), rep)
        if !haskey(JAM_DiFF_interpolators, key)
            if !haskey(dict_raw_DiFF, key)
                throw(KeyError(key))
            end
            build_JAM_DiFF_interpolator!(key)
        end
        check_JAM_DiFF_bounds(kind = kind, z = z, Mh = Mh, mu = mu, rep = rep)
        return JAM_DiFF_interpolators[key](Float32(z), Float32(Mh))
    end
end
