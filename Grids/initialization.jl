include("interpolation.jl")

@inline mu_key(mu::Real) = Float16(mu)

if @isdefined(dict_raw_DiFF) && @isdefined(selected_keys)
    global JAM_DiFF_interpolators = Dict{Any,Any}()

    zvals = Float64.(dict_raw_DiFF["z_array"])
    Mhvals = Float64.(dict_raw_DiFF["Mh_array"])
    flavors_by_kind = haskey(dict_raw_DiFF, "flavors_by_kind") ? dict_raw_DiFF["flavors_by_kind"] : nothing

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
        return JAM_DiFF_interpolators[key](Float32(z), Float32(Mh))
    end
end
