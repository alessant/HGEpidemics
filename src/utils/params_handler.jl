function initialize_params(config::Dict; int_max::Int = typemax(Int), df_ac::Union{DataFrame,Nothing} = nothing, printme::Bool = false)

    intervention_start = haskey(config, :intervention_start) ? config[:intervention_start] : int_max

    # PPM
    αₚ = haskey(config, :αₚ) ? config[:αₚ] : 0
    ppm_βd = haskey(config, :ppm_βd) ? config[:ppm_βd] : nothing
    ppm_βₑ = haskey(config, :ppm_βₑ) ? config[:ppm_βₑ] : nothing
    ppm_βᵢ = haskey(config, :ppm_βᵢ) ? config[:ppm_βᵢ] : nothing

    if haskey(config, :nodes_selection_strategy)
        nodes_selection_strategy =
            ismissing(config[:nodes_selection_strategy]) ? nothing : getfield(Main, Symbol(config[:nodes_selection_strategy]))
    else
        nodes_selection_strategy = nothing
    end

    # EM
    sanitize = haskey(config, :sanitize) ? config[:sanitize] : nothing
    sanification_intervals = haskey(config, :sanification_intervals) ? config[:sanification_intervals] : nothing

    # ISOLATION
    βisolation = haskey(config, :βisolation) ? config[:βisolation] : nothing

    #QUARANTINE
    βquarantine = haskey(config, :βquarantine) ? config[:βquarantine] : nothing

    # TRACING
    αᵢ = haskey(config, :αᵢ) ? config[:αᵢ] : 0
    βtracing = haskey(config, :βtracing) ? config[:βtracing] : 0

    if haskey(config, :nodes_tracing_strategy)
        nodes_tracing_strategy =
            ismissing(config[:nodes_tracing_strategy]) ? nothing : getfield(Main, Symbol(config[:nodes_tracing_strategy]))
    else
        nodes_tracing_strategy = nothing
    end

    # LOCKDOWN
    αₑ = haskey(config, :αₑ) ? config[:αₑ] : 0

    if haskey(config, :edges_selection_strategy)
        edges_selection_strategy =
            ismissing(config[:edges_selection_strategy]) ? nothing : getfield(Main, Symbol(config[:edges_selection_strategy]))
    else
        edges_selection_strategy = nothing
    end

    if haskey(config, :path)
        edges_selection_strategy_kwargs =
            ismissing(config[:path]) ? Dict{}() : Dict{Symbol, String}(:path => config[:path])
    else
        edges_selection_strategy_kwargs = Dict{}()
    end

    # AVOIDING CROWDING 
    df_avoiding_crowding = haskey(config, :avoiding_crowding) && config[:avoiding_crowding] ? df_ac : nothing

    thr_ac = haskey(config, :thr_ac) && !ismissing(config[:thr_ac]) ? config[:thr_ac] : nothing  

    #
    # params dict
    #
    params = OrderedDict{Symbol, Any}(
        :intervention_start => intervention_start,
        :αₚ => αₚ,
        :ppm_βd => ppm_βd, 
        :ppm_βₑ => ppm_βₑ, 
        :ppm_βᵢ => ppm_βᵢ,
        :nodes_selection_strategy => nodes_selection_strategy,
        :sanitize => sanitize,
        :sanification_intervals => sanification_intervals,
        :βisolation => βisolation,
        :βquarantine => βquarantine,
        :βtracing => βtracing,
        :αᵢ => αᵢ,
        :nodes_tracing_strategy => nodes_tracing_strategy,
        :αₑ => αₑ,
        :edges_selection_strategy => edges_selection_strategy,
        :edges_selection_strategy_kwargs => edges_selection_strategy_kwargs,
        :df_avoiding_crowding => df_avoiding_crowding,
        :thr_ac => thr_ac
    )

    printme && print_params(params)

    params
end


function print_params(params)
    _params = deepcopy(params)
    if !isnothing(params[:df_avoiding_crowding])
        _params[:df_avoiding_crowding] = "avoiding_crowding"
        println(_params)
    else
        println(params)
    end
end


