"""
    simulate(
            sim_type::SIS,
            df::DataFrame,
            intervals::Dict{Int, Pair{DateTime, DateTime}},
            user2vertex::Dict{String, Int},
            loc2he::Dict{String, Int},
            δ::Dates.Millisecond;
            Δ::Union{Int,TimePeriod,Nothing} = nothing,
            vstatus::Union{Array{Int, 1}, Nothing} = nothing,
            per_infected::Float64 = 0.2,
            c::Union{Int, Nothing} = 5,
            βd::Float64 = 0.2,
            βₑ::Float64 = 0.06,
            βᵢ::Float64 = 0.1,
            γₑ::Float64 = 0.06,
            γₐ::Float64 = 0.1,
            αᵥ::Union{Int, Float64} = 0.0,
            αₑ::Union{Int, Float64} = 0.0,
            lockdown::Bool = false,
            βₗ::Union{Int, Float64} = 0.0,
            intervention_start::Int = 0,
            nodes_imm_strategy::Union{Function, Nothing} = nothing,
            hes_imm_strategy::Union{Function, Nothing} = nothing,
            nodes_kwargs::Dict = Dict{}(),
            hes_kwargs::Dict = Dict{}(),
            niter::Int = 1,
            output_path::Union{AbstractString, Nothing} = nothing,
            kwargs...
    )

 Simulate a Susceptible-Infected-Susceptible diffusion model exploiting a
 Time-Varying Hypergraph. An immunization strategy may be applied either
 on vertices (i.e. people) and hyperedges (i.e. locations).

 **Note**
Sanite locations at regular time intervals.

 **Arguments**
 - `sim_type`, diffusion model to simulate;
 - `df`, DataFrame containing check-in data;
 - `intervals`, time intervals within which an indirect contact may happen;
 - `user2vertex`, mapping from agent ids to vertex ids;
 - `loc2he`, mapping from location ids to hyperedge ids;
 - `δ`, time within which a direct contact may happen (expressed in milliseconds);
 - `Δ`, time within which an indirect contact may happen (just for logging purposes);
 - `vstatus`, initial status (susceptible or infected) for each node. If it not
 given as input, it will be initialized in the simulation;
 - `per_infected`, percetage of initial infected agents;
 - `c`, factor bounding the probability to become infected when the number of
 contact increases;
 - `βd`, probability of becoming infected with a direct contact;
 - `βₑ`, probability that a location is infected by an agent;
 - `βᵢ`, probability of becoming infected with an indirect contact;
 - `γₑ`, probability that a location spontaneously recovers;
 - `γₐ`, probability that an agent spontaneously recovers;
 - `αᵥ`, percentage of agents to immunize;
 - `αₑ`, percentage of locations to santize;
 - `lockdown`, whether applying a lockdown policy;
 - `βₗ`, factor decreasing the probability that a sanitized location infects
 an agent;
 - `intervention_start`, iteration from which the immunization phase takes place;
 - `nodes_imm_strategy`, immunization strategy to apply to the agents;
 - `hes_imm_strategy`, immunization strategy to apply to the hyperedges;
 - `nodes_kwargs`, optional params for `nodes_imm_strategy`;
 - `hes_kwargs`, optional params for `hes_imm_strategy`;
 - `niter`, number of iteration the simulation is repeated;
 - `output_path`, path where logs are stored;
 - `kwargs`, other optional params.

"""
function simulate(
        sim_type::SIS_NPIs,
        df::DataFrame,
        #df_avoiding_crowding::DataFrame, #TODO modify in kwargs
        intervals::Dict{Int, Pair{DateTime, DateTime}},
        user2vertex::Dict{String, Int},
        loc2he::Dict{String, Int},
        δ::Dates.Millisecond;
        Δ::Union{Int,TimePeriod,Nothing} = nothing,
        vstatus::Union{Array{Int, 1}, Nothing} = nothing,
        per_infected::Float64 = 0.2,
        c::Union{Int, Nothing} = 5,
        βd::Union{Int, Float64} = 0.2,
        βₑ::Union{Int, Float64} = 0.06,
        βᵢ::Union{Int, Float64} = 0.1,
        γₑ::Union{Int, Float64} = 0.06,
        γₐ::Union{Int, Float64} = 0.1,
        #αᵥ::Union{Int, Float64} = 0.0, #TODO modify in kwargs
        #αₑ::Union{Int, Float64, Nothing} = 0.0, #TODO modify in kwargs
        #lockdown::Bool = false, #TODO modify in kwargs
        #βₗ::Union{Int, Float64} = 0.0, #TODO modify in kwargs
        #intervention_start::Int = 0, #TODO modify in kwargs
        #nodes_imm_strategy::Union{Function, Nothing} = nothing, #TODO modify in kwargs
        #hes_imm_strategy::Union{Function, Nothing} = nothing, #TODO modify in kwargs
        #nodes_kwargs::Dict = Dict{}(), #TODO modify in kwargs
        #hes_kwargs::Dict = Dict{}(), #TODO modify in kwargs
        niter::Int = 10,
        output_path::Union{AbstractString, Nothing} = nothing,
        print_me::Bool = true,
        store_me::Bool = true,
        kwargs...
)

    #########################
    # Check NPIs params
    ########################
    df_avoiding_crowding = haskey(kwargs, :df_avoiding_crowding) ? kwargs[:df_avoiding_crowding] : nothing
    avoiding_crowding = isnothing(df_avoiding_crowding) ? false : true
    thr_ac = haskey(kwargs, :thr_ac) ? kwargs[:thr_ac] : nothing

    βisolation = haskey(kwargs, :βisolation) ? kwargs[:βisolation] : nothing
    βquarantine = haskey(kwargs, :βquarantine) ? kwargs[:βquarantine] : nothing
    β_quarantine_loc = haskey(kwargs, :β_quarantine_loc) ? kwargs[:β_quarantine_loc] : nothing

    αₚ = haskey(kwargs, :αₚ) ? kwargs[:αₚ] : 0
    ppm = αₚ > 0 ? true : false
    ppm_βₑ = haskey(kwargs, :ppm_βₑ) ? kwargs[:ppm_βₑ] : nothing
    ppm_βd = haskey(kwargs, :ppm_βd) ? kwargs[:ppm_βd] : nothing
    ppm_βᵢ = haskey(kwargs, :ppm_βᵢ) ? kwargs[:ppm_βᵢ] : nothing

    αᵥ = haskey(kwargs, :αᵥ) ? kwargs[:αᵥ] : 0
    nodes_immunization = αᵥ > 0 ? true : false

    αₑ = haskey(kwargs, :αₑ) ? kwargs[:αₑ] : 0
    lockdown = αₑ > 0 ? true : false

    αᵢ = haskey(kwargs, :αᵢ) ? kwargs[:αᵢ] : 0
    tracing = αᵢ > 0 ? true : false
    βtracing = haskey(kwargs, :βtracing) ? kwargs[:βtracing] : 0

    βₗ = haskey(kwargs, :βₗ) ? kwargs[:βₗ] : 0
    sanitize = haskey(kwargs, :sanitize) ? kwargs[:sanitize] : nothing
    sanification_intervals = haskey(kwargs, :sanification_intervals) ? kwargs[:sanification_intervals] : Array{Int, 1}()
    # println(sanification_intervals)

    intervention_start = haskey(kwargs, :intervention_start) ? kwargs[:intervention_start] : length(intervals) + 1

    nodes_selection_strategy = haskey(kwargs, :nodes_selection_strategy) ? kwargs[:nodes_selection_strategy] : nothing
    nodes_tracing_strategy = haskey(kwargs, :nodes_tracing_strategy) ? kwargs[:nodes_tracing_strategy] : nothing
    edges_selection_strategy = haskey(kwargs, :edges_selection_strategy) ? kwargs[:edges_selection_strategy] : nothing
    
    nodes_selection_strategy_kwargs = haskey(kwargs, :nodes_selection_strategy_kwargs) ? kwargs[:nodes_selection_strategy_kwargs] : Dict{}()
    nodes_tracing_strategy_kwargs = haskey(kwargs, :nodes_tracing_strategy_kwargs) ? kwargs[:nodes_tracing_strategy_kwargs] : Dict{}()
    edges_selection_strategy_kwargs = haskey(kwargs, :edges_selection_strategy_kwargs) ? kwargs[:edges_selection_strategy_kwargs] : Dict{}()

    #########################
    # Init logs
    ########################
    if store_me 
        if isnothing(output_path)
            if !isdir("results/")
                Filesystem.mkdir("results/")
            end
            output_path = "results/$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv"
        end

        header = string(
            "sim_step,Δ,δ,c,per_infected,βd,βₑ,βᵢ,γₑ,γₐ,",
            "avg_he_size,avg_degree,avg_direct_contacts,new_agents,",
            "moved_agents,perc_infected_agents,perc_infected_locations,",
            "avoiding_crowding,αᵥ,αₑ,lockdown,βₗ,intervention_start",
            "nodes_selection_strategy,edges_selection_strategy\n"
            )

        init_log(output_path, header)
    end

    #########################
    # Init simulation params
    ########################

    # iter -> percentage of infected agents per simulation step
    to_return = Dict{Int, Array{Float64, 1}}()
    # iter -> absolute number of infected agents per simulation step
    to_return_abs = Dict{Int, Array{Int, 1}}()

    # new infected x day due to direct contact
    new_infected_direct_to_return = Dict{Int, Array{Int, 1}}()
    # new infected x day due to indirect contact
    new_infected_indirect_to_return = Dict{Int, Array{Int, 1}}()
    # new infected x day due to either direct or indirect contact (sum of the previous two)
    new_infected_both_to_return = Dict{Int, Array{Int, 1}}()

    # number of isolated agents, if isolation measures are simulated
    isolated_to_return = Dict{Int, Array{Float64, 1}}()
    # number of quarantines agents, if quarantine measures are simulated
    quarantined_to_return = Dict{Int, Array{Float64, 1}}()
    # number of sick agents free to move (i.e. sick but not isolated)
    infected_not_isolated_to_return = Dict{Int, Array{Float64, 1}}()
    # number of sick agents free to move (i.e. sick but not quarantined)
    infected_not_quarantined_to_return = Dict{Int, Array{Float64, 1}}()
    # number of susceptible agents
    susceptible_to_return = Dict{Int, Array{Float64, 1}}()


    # agents met
    agents_met_to_return = Dict{Int, Array{Int, 1}}()
    # visited locations
    visited_locations_to_return = Dict{Int, Array{Int, 1}}()

    # evaluation of an initial vector of infected agents
    # if it is not given as input
    if isnothing(vstatus)
        vstatus = fill(1, length(user2vertex))
        vrand = rand(Float64, length(user2vertex))
        for i=1:length(user2vertex)
            if per_infected  <= vrand[i]
                vstatus[i] = 0
            end
        end
    end

    #########################
    # Start simulation 
    ########################

    # for randomization purposes
    for iter=1:niter
        print_me && println("Iter $(iter)")

        h = nothing
        added, moved = 0, 0

        # percentage of infected per simulation step
        per_infected_sim = Array{Float64, 1}()
        n_infected = Array{Int, 1}()

        new_infected_direct = Array{Int, 1}()
        new_infected_indirect = Array{Int, 1}()
        new_infected_both = Array{Int, 1}()

        # number of isolated and quarantined 
        isolated_users = Array{Float64, 1}()
        quarantined_users = Array{Float64, 1}()
        infected_not_isolated = Array{Float64, 1}()
        infected_not_quarantined = Array{Float64, 1}()
        susceptible = Array{Float64, 1}()

        # agents seen 
        agents_met = fill(0, length(user2vertex), length(user2vertex))
        # location visited 
        location_visited = fill(0, length(user2vertex), length(loc2he))

        # store which agents are present in the given timeframe
        agentsepoc = zeros(Int, length(user2vertex))

        # Initialize the status vector of the nodes
        _vstatus = copy(vstatus)
        # Initially, all location are safe
        hestatus = zeros(Int, length(loc2he))

        # Storing the new status of each vertex and hyperedge
        vnextstatus = copy(vstatus)
        henextstatus = copy(hestatus)

        n_users = length(_vstatus)

        push!(per_infected_sim, sum(_vstatus) / n_users)
        push!(n_infected, sum(_vstatus))

        push!(isolated_users, sum(isolated_users) / length(_vstatus))
        push!(quarantined_users, sum(quarantined_users) / length(_vstatus))
        push!(susceptible, (length(_vstatus) - sum(_vstatus)) / length(_vstatus))

        ################
        # IMMUNIZATION
        # the agent can no longer contract the infection
        # the same happens for a location
        # however, they can still spread it
        ################
        istatus = fill(0, length(user2vertex))
        ihestatus = fill(0, length(loc2he))

        nextistatus = copy(istatus)
        nextihestatus = copy(ihestatus)

        ######################à
        # PPMs
        # store which agent is using PPMs
        #####################
        ppm_agents = fill(0, length(user2vertex))

        ######################à
        # SDM 
        # - isolation
        # - quarantine
        #####################
        isolation = fill(0, length(user2vertex))
        quarantine = fill(0, length(user2vertex))
        
        next_isolation = fill(0, length(user2vertex))
        next_quarantine = fill(0, length(user2vertex))

        push!(infected_not_isolated, (sum(_vstatus) - sum(isolation)) / length(_vstatus))
        push!(infected_not_quarantined, (sum(_vstatus) - sum(quarantine)) / length(_vstatus))

        ################
        # TRACING
        ################
        to_trace = Dict{Int, Set{Int}}()
        users_to_trace = Set{Int}()

        ################
        # SIMULATION
        ################
        for t=1:length(intervals)
            checkins_df = avoiding_crowding && t >= intervention_start ? df_avoiding_crowding : df
            thr = !isnothing(thr_ac) && t >= intervention_start ? thr_ac : nothing

            h, added, moved = generatehg!(
                    h,
                    checkins_df,
                    get(intervals, t, 0).first,
                    get(intervals, t, 0).second,
                    user2vertex,
                    loc2he,
                    agentsepoc,
                    t;
                    thr = thr
                )

            isnothing(h) && continue

            new_infected_dir = 0
            new_infected_ind = 0

            if !isnothing(thr_ac) && t > intervention_start
                for he in 1:nhe(h)
                    if length(getvertices(h, he)) > thr_ac
                        println("Over thr ", length(getvertices(h, he)))
                    end
                end
            end

            # Estimation of the parameter c
            # based on the distribution
            # of the hyperedge size
            if isnothing(c)
                dist = Array{Int, 1}()

                for he=1:nhe(h)
                    push!(dist, length(getvertices(h, he)))
                end

                c = median(dist)
                println(t, " -- ", c)
            end


            #################################
            # Evaluating some stats for
            # the current hg
            #################################

            # hyperedges average size
            avg_he_size = .0
            for he=1:nhe(h)
                avg_he_size += length(getvertices(h, he))
            end
            avg_he_size /= nhe(h)

            # nodes average degree
            avg_degree = .0
            for v=1:nhv(h)
                avg_degree += length(gethyperedges(h, v))
            end
            avg_degree /= nhv(h)

            # number of infected locations with
            # at least two agents
            infected_locations = 0
            for he=1:nhe(h)
                if hestatus[he] == 1 && length(getvertices(h, he)) > 1
                    infected_locations += 1
                end
            end


            ########################
            # IMMUNIZATION
            # immunizing a node does not have effect on its S/I status:
            # a node that has been immunized cannot get sick,
            # but it can still spread the contagion
            ########################

            # start the immunization process
            # in a given timeframe
            if t == intervention_start
                # apply the given selection strategy
                # to choose which αₚ nodes will use PPM
                if ppm
                    using_ppms = nodes_selection_strategy(h, αₚ; nodes_selection_strategy_kwargs...)
                    map(v -> ppm_agents[v] = 1, using_ppms)
                end

                # apply the given selection strategy
                # to choose which αᵥ nodes will be immunized
                if nodes_immunization
                    to_immunize = nodes_selection_strategy(h, αᵥ; nodes_selection_strategy_kwargs...)
                    map(v -> nextistatus[v] = 1, to_immunize)
                end

                # apply the given immunization strategy
                # over αₑ hyperedges
                # LOCKDOWN
                if lockdown
                    to_close = edges_selection_strategy(dual(h), αₑ; edges_selection_strategy_kwargs...)
                    map(he -> nextihestatus[he] = 1, to_close)
                end

                if tracing
                    u2trace = nodes_tracing_strategy(h, αᵢ; nodes_tracing_strategy_kwargs...)
                    map(user -> to_trace[user] = Set{Int}(), u2trace)
                    users_to_trace = keys(to_trace)
                end
            end

            ########################
            # ISOLATION
            # agents may decide to isolate theirself
            # if their are sick and not immunized or already quarantined_percentage
            # according to the probability defined by βisolation
            ########################
            if !isnothing(βisolation) && t >= intervention_start
                for v=1:nhv(h)
                    if agentsepoc[v] == 1
                        if _vstatus[v] == 1 && istatus[v] == 0 && quarantine[v] == 0 && rand() < 1 -  ℯ ^ - (kwargs[:βisolation])
                            next_isolation[v] = 1
                        end
                    end
                end
            end

            ########################
            # TRACING
            # If an agent is not in quarantined yet,
            # it may enter the quarantine with a probability
            # proportional to the number of infected contacts.
            ########################
            if tracing && t >= intervention_start
                for v in users_to_trace
                    if quarantine[v] == 0
                        i = 0
                        for u in to_trace[v]
                            i += _vstatus[u]
                        end
                        if rand() < 1 - ℯ ^ - (βtracing * i)
                            quarantine[v] = 1
                        end
                    end
                end
            end

            ########################
            # DIFFUSION ALGORITHM
            ########################

            #
            # PHASE 1 - Agent-to-Environment
            #
            for he=1:nhe(h)
                # If the location is immunized
                # and the lockdown is active,
                # it cannot spread the infection anymore
                #! with this implementation
                #! the value of lockdown depends on αₑ
                ihestatus[he] == 1 && lockdown && continue 

                # If the location has at least two agents
                # and it is not infected, it may become contamined.
                if length(getvertices(h, he)) > 1 && hestatus[he] == 0
                    prob = nothing

                    # if there are agents that are using ppms
                    if ppm
                        i, ppm_i = infected_ppm(h, he, _vstatus, ppm_agents; istatus = istatus, quarantine = quarantine, isolation = isolation)
                        prob = 1 - ℯ ^ - (βₑ * f(i, c) + ppm_βₑ * f(ppm_i, c))
                    else
                        #println("argh 1")
                        i = infected(h, he, _vstatus; istatus = istatus, quarantine = quarantine, isolation = isolation)
                        prob = 1 - ℯ ^ - (βₑ * f(i, c))
                    end

                    # i, ppm_i = infected_ppm(h, he, _vstatus, istatus, ppm_agents)
                    # prob = 1 - ℯ ^ - (βₑ * f(i, c) + kwargs[:ppm_βₑ] * f(ppm_i, c))

                    # i = infected(h, he, _vstatus, istatus; quarantine = quarantine, isolation = isolation)
                    # prob2 = 1 - ℯ ^ - (βₑ * f(i, c))   
                    
                    # println("PPM prob = $prob --- NO PPM prob = $prob2")

                    if rand() < prob
                        # the location is contaminated
                        henextstatus[he] = 1
                    end
                elseif rand() <  1 - ℯ ^ - γₑ
                    # the location has been sanitized
                    henextstatus[he] = 0
                end
            end

            #
            # PHASE 2 - Agent-to-Agent
            #
            avg_direct_contacts = 0
            for v=1:nhv(h)
                # a node that has been immunized or it is isolated or quarantined
                # cannot become infected
                (istatus[v] == 1 || isolation[v] == 1 || quarantine[v] == 1) && continue

                # if the agent is present in the current timeframe
                if agentsepoc[v] == 1
                    i = 0
                    for he in gethyperedges(h, v)

                        # if the lockdown is active,
                        # a direct contact in that location cannot happen
                        #! same as phase 1
                        ihestatus[he.first] == 1 && lockdown && continue

                        for u in getvertices(h, he.first)
                            if v != u.first
                                # if u and v have been together in the same place
                                # in a time interval less than δ
                                # then it counts ad a direct contact
                                if abs(h[v, he.first] - h[u.first, he.first]) <= δ.value

                                    # v and u have met 
                                    # after the application of an intervention
                                    if t > intervention_start
                                        agents_met[v, u.first] = 1
                                    end

                                    # if v is healthy and
                                    # u is not quarantines/isolated/immunized
                                    if _vstatus[v] == 0 && isolation[u.first] == 0 && quarantine[u.first] == 0 && istatus[u.first] == 0
                                        i += _vstatus[u.first]

                                        if tracing #haskey(kwargs, :tracing) && !isnothing(kwargs[:tracing])
                                            if v in users_to_trace && u.first in users_to_trace
                                                push!(
                                                    get!(to_trace, v, Set{Int}()),
                                                    u.first
                                                )
        
                                                push!(
                                                    get!(to_trace, u.first, Set{Int}()),
                                                    v
                                                )
                                            end
                                        end
                                    end
                                    avg_direct_contacts += 1
                                end
                            end
                        end
                    end
                    #if an agent is using PPMs, then
                    #the probability to be infected from
                    #another agent is less
                    if ppm_agents[v] == 1 
                        if _vstatus[v] == 0 && rand() < 1 - ℯ ^ - (ppm_βd * i)
                            #println("$v using ppm")
                            vnextstatus[v] = 1
                            new_infected_dir += 1
                        end
                    # an agent becomes infected according to
                    # the following probability
                    else
                        if _vstatus[v] == 0 && rand() < 1 - ℯ ^ - (βd * i)
                            #println("argh 2 -- $v -- using ppm $(ppm_agents[v])")
                            vnextstatus[v] = 1
                            new_infected_dir += 1
                        end
                    end

                    # go to quarantine
                    if !isnothing(βquarantine) && t >= intervention_start
                        if rand() < 1 - ℯ ^ - (βquarantine * i)
                            next_quarantine[v] = 1
                        end
                    end
                end
            end

            avg_direct_contacts \= sum(agentsepoc)


            #
            # PHASE 3 - Environment-to-Agent
            #
            for v=1:nhv(h)

                # if the agent is present in the current timeframe
                if agentsepoc[v] == 1

                    for he in gethyperedges(h, v)
                        #! same as phase 1
                        # if the location is closed
                        # I cannot checking in 
                        ihestatus[he.first] == 1 && lockdown && continue
                        # if the agent is in isolation
                        # or in quarantine, it cannot checkin
                        (isolation[v] == 1 || quarantine[v] == 1) && continue

                        # store that an agent has been in that location
                        # after the protective measure
                        if t > intervention_start
                            location_visited[v, he.first] = 1
                        end
                    end

                    # if the agent is healthy and
                    # it is not immunized/quarantined/isolated,
                    # it may become infected
                    if _vstatus[v] == 0 && istatus[v] == 0 && isolation[v] == 0 && quarantine[v] == 0
                        i = 0
                        i_immunized = 0

                        for he in gethyperedges(h, v)

                            # if the lockdown is active,
                            # an indirect contact in that location
                            # cannot take place
                            #! same as phase 1
                            ihestatus[he.first] == 1 && lockdown && continue

                            # store that an agent has been in that location
                            # if t >= intervention_start
                            #     location_visited[v, he.first] = 1
                            # end

                            if length(getvertices(h, he.first)) > 1
                                if ihestatus[he.first] == 0
                                    i += hestatus[he.first]
                                else
                                    i_immunized += hestatus[he.first]
                                end
                            end
                        end

                        #if an agent is using a PPM, then
                        #the probability to be infected from
                        #an environment is less
                        if ppm_agents[v] == 1 
                            if rand() < 1 - ℯ ^ -( (ppm_βᵢ * i) + ( (ppm_βᵢ * (1 - βₗ)) * i_immunized) )
                                vnextstatus[v] = 1
                                new_infected_ind += 1
                            end
                        else
                            if rand() < 1 - ℯ ^ -( (βᵢ * i) + ( (βᵢ * (1 - βₗ)) * i_immunized) ) #1 - ℯ^-(βᵢ * f(i, c))
                                #println("argh 3 -- $v -- using ppm $(ppm_agents[v])")
                                vnextstatus[v] = 1
                                new_infected_ind += 1
                            end
                        end

                        # go to quarantine
                        # according to the visited locations
                        if !isnothing(β_quarantine_loc)
                            if rand() < 1 - ℯ ^ - (β_quarantine_loc * i)
                                next_quarantine[v] = 1
                            end
                        end

                    elseif rand() < 1 - ℯ ^ - γₐ
                        # the agent spontaneously returns healthy
                        vnextstatus[v] = 0
                        next_isolation[v] = 0
                        next_quarantine[v] = 0
                    end
                end
            end

            #TODO implement log

            push!(per_infected_sim, sum(_vstatus) / n_users)
            push!(n_infected, sum(_vstatus))

            push!(new_infected_direct, new_infected_dir)
            push!(new_infected_indirect, new_infected_ind)
            push!(new_infected_both, (new_infected_dir+new_infected_ind))

            push!(isolated_users, sum(isolation)  / length(_vstatus))
            push!(quarantined_users, sum(quarantine)  / length(_vstatus))
            push!(infected_not_isolated, (sum(_vstatus) - sum(isolation)) / length(_vstatus))
            push!(infected_not_quarantined, (sum(_vstatus) - sum(quarantine)) / length(_vstatus))
            push!(susceptible, (length(_vstatus) - sum(_vstatus)) / length(_vstatus))
            
            # update
            _vstatus = copy(vnextstatus)
            istatus = copy(nextistatus)

            hestatus = copy(henextstatus)
            ihestatus = copy(nextihestatus)

            isolation = copy(next_isolation)
            quarantine = copy(next_quarantine)

            # Δ=4
            # At the end of each day,
            # all locations are sanitized

            # Δ = 12,
            # All locations are sanitized
            # at midday
            if !isnothing(sanitize) && sanitize && t >= intervention_start
                if t % (convert(Int, 24 / Δ)) in sanification_intervals #[3, 4, 5] BLE 
                    #println("Sanitizing $(sum(hestatus)) @ interval $(t) -- $(get(intervals, t, 0).first) - $(get(intervals, t, 0).second)")
                    hestatus = zeros(Int, length(loc2he))
                end
            end
        end

        push!(to_return, iter=>per_infected_sim)
        push!(to_return_abs, iter=>n_infected)

        push!(new_infected_direct_to_return, iter=>new_infected_direct)
        push!(new_infected_indirect_to_return, iter=>new_infected_indirect)
        push!(new_infected_both_to_return, iter=>new_infected_both)

        push!(isolated_to_return, iter => isolated_users)
        push!(quarantined_to_return, iter => quarantined_users)
        push!(infected_not_isolated_to_return, iter => infected_not_isolated)
        push!(infected_not_quarantined_to_return, iter => infected_not_quarantined)
        push!(susceptible_to_return, iter => susceptible)

        push!(agents_met_to_return, iter => vec(sum(agents_met, dims=1)))
        push!(visited_locations_to_return, iter => vec(sum(location_visited, dims=2)))
    end
    

    Dict{Symbol, Any}(
        :infected_percentage => to_return,
        :n_infected => to_return_abs,
        :new_infected_direct => new_infected_direct_to_return,
        :new_infected_indirect => new_infected_indirect_to_return,
        :new_infected_both => new_infected_both_to_return,
        :isolated_percentage => isolated_to_return,
        :quarantined_percentage => quarantined_to_return,
        :infected_not_isolated_percentage => infected_not_isolated_to_return,
        :infected_not_quarantined_percentage => infected_not_quarantined_to_return,
        :susceptible_percentage => susceptible_to_return,
        :agents_met => agents_met_to_return,
        :location_visited => visited_locations_to_return
    )
end
