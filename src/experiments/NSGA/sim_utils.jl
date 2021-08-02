function init_simulation(data_params)
    input_data = JSON.parse((open(data_params, "r")))
    
    fdata_params = input_data["data_params"]
    fparams = input_data["sim_params"]
    
    jtable = jsontable(read(open(fparams, "r")))
    paramsdf = DataFrame(jtable)
    
    #########################
    # Generating model data
    ########################
    data_params = JSON3.read(read(open(fdata_params, "r")))
    header = [Symbol(col) for col in data_params.header]

    df, intervals, user2vertex, loc2he =
        generate_model_data(
            data_params.dataset,
            header,
            Symbol(data_params.userid),
            Symbol(data_params.venueid),
            Symbol(data_params.UTCtime),
            data_params.dateformat;
            datarow = data_params.datarow,
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf.Δ[1])),
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf.δ[1])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )

    df_ac, _, _, _ =
        generate_model_data(
            data_params.dataset_ac,
            header,
            Symbol(data_params.userid),
            Symbol(data_params.venueid),
            Symbol(data_params.UTCtime),
            data_params.dateformat_ac;
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf.Δ[1])),
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf.δ[1])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date),
            datarow = data_params.datarow_ac
        )
    
    #########################
    # Initialization of infected nodes
    ########################
 
    users = keys(user2vertex)
    
    # vstatus = fill(1, length(users))
    # vrand = rand(Float64, (1, length(users)))

    # for i=1:length(users)
    #     if infected_percentage <= vrand[i]
    #         vstatus[i] = 0
    #     end
    # end

    vstatus = fill(0, length(users))
    vstatus[1] = 1

    # check lockdown strategy
    edges_selection_strategy = 
        "edges_selection_strategy" in names(paramsdf) ? getfield(Main, Symbol(paramsdf.edges_selection_strategy[1])) : uniform

    edges_selection_strategy_kwargs = 
        "path" in names(paramsdf) ? Dict{Symbol, String}(:path => paramsdf.path[1]) : Dict{}()

    (
        df = df,
        df_ac = df_ac,
        intervals = intervals,
        user2vertex = user2vertex,
        loc2he = loc2he,
        vstatus = vstatus,
        Δ = paramsdf.Δ[1],
        δ = paramsdf.δ[1],
        βd = paramsdf.βd[1],
        βᵢ = paramsdf.βᵢ[1],
        βₑ = paramsdf.βₑ[1],
        γₑ = paramsdf.γₑ[1],
        γₐ = paramsdf.γₐ[1],
        niter = paramsdf.niter[1],
        intervention_start = paramsdf.intervention_start[1],
        ppm_βd = paramsdf.ppm_βd[1],
        ppm_βₑ = paramsdf.ppm_βₑ[1],
        ppm_βᵢ = paramsdf.ppm_βᵢ[1],
        sanification_intervals = paramsdf.sanification_intervals[1],
        βtracing = paramsdf.βtracing[1],
        edges_selection_strategy = edges_selection_strategy,
        edges_selection_strategy_kwargs = edges_selection_strategy_kwargs
    )
end


function run_simulation(config::Vector{Float64}, simulation_data::NamedTuple; printme::Bool=false)
    params = parse_conf(config)

    npi_paramas = Dict{}(
            :intervention_start => simulation_data.intervention_start,
            :αₚ => params[:αₚ],
            :ppm_βd => simulation_data.ppm_βd, 
            :ppm_βₑ => simulation_data.ppm_βₑ, 
            :ppm_βᵢ => simulation_data.ppm_βᵢ,
            :nodes_selection_strategy => uniform,
            :sanitize => params[:sanitize],
            :sanification_intervals => simulation_data.sanification_intervals,
            :βisolation => params[:βisolation],
            :βquarantine => params[:βquarantine],
            :αₑ => params[:αₑ],
            :βtracing => simulation_data.βtracing,
            :αᵢ => params[:αᵢ],
            :nodes_tracing_strategy => uniform,
            :edges_selection_strategy => simulation_data.edges_selection_strategy,
            :edges_selection_strategy_kwargs => simulation_data.edges_selection_strategy_kwargs,
            :df_avoiding_crowding => params[:avoiding_crowding] ? simulation_data.df_ac : nothing,
        )

    printme && println(npi_paramas)

    results =
        simulate(
            SIS_NPIs(),
            simulation_data.df,
            simulation_data.intervals,
            simulation_data.user2vertex,
            simulation_data.loc2he,
            convert(Dates.Millisecond, Dates.Minute(simulation_data.δ));
            Δ = simulation_data.Δ,
            vstatus = simulation_data.vstatus,
            βd = simulation_data.βd,
            βᵢ = simulation_data.βᵢ,
            βₑ = simulation_data.βₑ,
            γₑ = simulation_data.γₑ,
            γₐ = simulation_data.γₐ,
            niter = simulation_data.niter,
            store_me = false,
            print_me = false,
            npi_paramas...
        )

    # get the average over each simulation run
    infected_distribution = mean(collect(values(results[:infected_percentage])))
    susceptible_distribution = mean(collect(values(results[:susceptible_percentage])))
    isolated_distribution = mean(collect(values(results[:isolated_percentage])))
    quarantined_distribution = mean(collect(values(results[:quarantined_percentage])))

    infected_not_isolated_distribution = mean(collect(values(results[:infected_not_isolated_percentage])))
    infected_not_quarantined_distribution = mean(collect(values(results[:infected_not_quarantined_percentage])))

    agents_met = mean(collect(values(results[:agents_met])))
    location_visited = mean(collect(values(results[:location_visited])))

    (
        infected_distribution = infected_distribution, 
        susceptible_distribution = susceptible_distribution,
        isolated_distribution = isolated_distribution,
        quarantined_distribution = quarantined_distribution,
        infected_not_isolated_distribution = infected_not_isolated_distribution,
        infected_not_quarantined_distribution = infected_not_quarantined_distribution,   
        agents_met = agents_met,
        location_visited = location_visited,
    )
end


function parse_conf(config)
    params = Dict{Symbol, Union{Float64, Bool}}()

    params[:αₚ] = config[1]
    params[:sanitize] = config[2] >= 0.5 ? true : false
    params[:βisolation] = config[3]
    params[:βquarantine] = config[4]
    params[:αₑ] = config[5]
    params[:αᵢ] = config[6]
    params[:avoiding_crowding] = config[7] >= 0.5 ? true : false

    params
end
