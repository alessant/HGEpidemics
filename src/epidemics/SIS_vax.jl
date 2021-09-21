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
        sim_type::SIS_vax,
        df::DataFrame,
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
        βₗ::Union{Int, Float64} = 0.0,
        γₑ::Union{Int, Float64} = 0.06,
        γₐ::Union{Int, Float64} = 0.1,
        niter::Int = 10,
        output_path::Union{AbstractString, Nothing} = nothing,
        print_me::Bool = true,
        store_me::Bool = true,
        kwargs...
)

    #########################
    # Check immunization params
    ########################
    αᵥ = haskey(kwargs, :αᵥ) ? kwargs[:αᵥ] : 0
    nodes_immunization = αᵥ > 0 ? true : false

    αₑ = haskey(kwargs, :αₑ) ? kwargs[:αₑ] : 0
    lockdown = αₑ > 0 ? true : false

    intervention_start = haskey(kwargs, :intervention_start) ? kwargs[:intervention_start] : length(intervals) + 1

    # immunization params - nodes
    nodes_immunization_strategy = haskey(kwargs, :nodes_immunization_strategy) ? kwargs[:nodes_immunization_strategy] : nothing
    nodes_immunization_strategy_kwargs = haskey(kwargs, :nodes_immunization_strategy_kwargs) ? kwargs[:nodes_immunization_strategy_kwargs] : Dict{}()

    # immunization params - edges
    edges_immunization_strategy = haskey(kwargs, :edges_immunization_strategy) ? kwargs[:edges_immunization_strategy] : nothing
    edges_immunization_strategy_kwargs = haskey(kwargs, :edges_immunization_strategy_kwargs) ? kwargs[:edges_immunization_strategy_kwargs] : Dict{}()
    edges_immunization = isnothing(edges_immunization_strategy) ? false : true


    #########################
    # Init logs
    ########################
    if store_me 
        if isnothing(output_path)
            if !isdir("results/")SIS_vax
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

    # number of susceptible agents
    susceptible_to_return = Dict{Int, Array{Float64, 1}}()

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

        # number susceptible
        susceptible = Array{Float64, 1}()

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

        push!(susceptible, (length(_vstatus) - sum(_vstatus)) / length(_vstatus))

        ################
        # IMMUNIZATION
        # the agent can no longer contract or spread the infection
        # the same happens for a location
        ################
        istatus = fill(0, length(user2vertex))
        ihestatus = fill(0, length(loc2he))

        nextistatus = copy(istatus)
        nextihestatus = copy(ihestatus)

        ################
        # SIMULATION
        ################
        for t=1:length(intervals)

            h, added, moved = generatehg!(
                    h,
                    df,
                    get(intervals, t, 0).first,
                    get(intervals, t, 0).second,
                    user2vertex,
                    loc2he,
                    agentsepoc,
                    t
                )

            isnothing(h) && continue

            new_infected_dir = 0
            new_infected_ind = 0

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

            if t < intervention_start 
                #push!(per_infected_sim, sum(_vstatus) / (length(_vstatus) - sum(istatus)))
                push!(per_infected_sim, 0)

                push!(new_infected_direct, 0)
                push!(new_infected_indirect, 0)
                push!(new_infected_both, 0)

                push!(susceptible, 1)
                
                continue
            end

            # start the immunization process
            # in a given timeframe
            if t == intervention_start
                # apply the given immunization strategy
                # to choose which αᵥ nodes will be immunized
                if nodes_immunization

                    nodes_immunization_strategy_kwargs[:mindate] = get(intervals, kwargs[:start_slot], 0).first
                    nodes_immunization_strategy_kwargs[:maxdate] = get(intervals, kwargs[:end_slot], 0).second 

                    nodes_immunization_strategy_kwargs[:df] = df
                    nodes_immunization_strategy_kwargs[:intervals] = intervals
                    nodes_immunization_strategy_kwargs[:user2vertex] = user2vertex
                    nodes_immunization_strategy_kwargs[:loc2he] = loc2he
                    nodes_immunization_strategy_kwargs[:δ] = δ

                    nodes_immunization_strategy_kwargs[:imm_start] = intervention_start
                    nodes_immunization_strategy_kwargs[:start_slot] = kwargs[:start_slot]
                    nodes_immunization_strategy_kwargs[:end_slot] = kwargs[:end_slot]

                    to_immunize = nodes_immunization_strategy(h, αᵥ; nodes_immunization_strategy_kwargs...)
                    map(v -> nextistatus[v] = 1, to_immunize)
                end

                # apply the given immunization strategy
                # to choose which αₑ hyperedges will be immunized
                if edges_immunization

                    edges_immunization_strategy_kwargs[:mindate] = get(intervals, kwargs[:start_slot], 0).first
                    edges_immunization_strategy_kwargs[:maxdate] = get(intervals, kwargs[:end_slot], 0).second 

                    edges_immunization_strategy_kwargs[:df] = df
                    edges_immunization_strategy_kwargs[:intervals] = intervals
                    edges_immunization_strategy_kwargs[:user2vertex] = user2vertex
                    edges_immunization_strategy_kwargs[:loc2he] = loc2he
                    edges_immunization_strategy_kwargs[:δ] = δ

                    edges_immunization_strategy_kwargs[:imm_start] = intervention_start
                    edges_immunization_strategy_kwargs[:start_slot] = kwargs[:start_slot]
                    edges_immunization_strategy_kwargs[:end_slot] = kwargs[:end_slot]

                    to_immunize = edges_immunization_strategy(dual(h), αₑ; edges_immunization_strategy_kwargs...)
                    map(he -> nextihestatus[he] = 1, to_immunize)
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
                ihestatus[he] == 1 && continue 

                # If the location has at least two agents
                # and it is not infected, it may become contamined.
                if length(getvertices(h, he)) > 1 && hestatus[he] == 0
                    i = infected(h, he, _vstatus; istatus = istatus)
                    if rand() < 1 - ℯ ^ - (βₑ * f(i, c))
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
                istatus[v] == 1  && continue

                # if the agent is present in the current timeframe
                if agentsepoc[v] == 1
                    i = 0
                    for he in gethyperedges(h, v)

                        for u in getvertices(h, he.first)
                            if v != u.first
                                # if u and v have been together in the same place
                                # in a time interval less than δ
                                # then it counts ad a direct contact
                                if abs(h[v, he.first] - h[u.first, he.first]) <= δ.value

                                    # if v is healthy and
                                    # u is not immunized
                                    if _vstatus[v] == 0 && istatus[u.first] == 0
                                        i += _vstatus[u.first]
                                    end
                                    avg_direct_contacts += 1
                                end
                            end
                        end
                    end
                    # an agent becomes infected according to
                    # the following probability
                    if _vstatus[v] == 0 && rand() < 1 - ℯ ^ - (βd * i)
                        #println("argh 2 -- $v -- using ppm $(ppm_agents[v])")
                        vnextstatus[v] = 1
                        new_infected_dir += 1
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
                        ihestatus[he.first] == 1 && continue
                    end

                    # if the agent is healthy and
                    # it is not immunized/quarantined/isolated,
                    # it may become infected
                    if _vstatus[v] == 0 && istatus[v] == 0 
                        i = 0
                        i_immunized = 0

                        for he in gethyperedges(h, v)

                            # if the lockdown is active,
                            # an indirect contact in that location
                            # cannot take place
                            #! same as phase 1
                            ihestatus[he.first] == 1 && continue

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


                        if rand() < 1 - ℯ ^ -( (βᵢ * i) + ( (βᵢ * (1 - βₗ)) * i_immunized) ) #1 - ℯ^-(βᵢ * f(i, c))
                            vnextstatus[v] = 1
                            new_infected_ind += 1
                        end

                    elseif rand() < 1 - ℯ ^ - γₐ
                        # the agent spontaneously returns healthy
                        vnextstatus[v] = 0
                    end
                end
            end

            #TODO implement log

            push!(per_infected_sim, sum(_vstatus) / n_users)
            push!(n_infected, sum(_vstatus))

            push!(new_infected_direct, new_infected_dir)
            push!(new_infected_indirect, new_infected_ind)
            push!(new_infected_both, (new_infected_dir+new_infected_ind))

            push!(susceptible, (length(_vstatus) - sum(_vstatus)) / length(_vstatus))
            
            # update
            _vstatus = copy(vnextstatus)
            istatus = copy(nextistatus)

            hestatus = copy(henextstatus)
            ihestatus = copy(nextihestatus)
        end

        push!(to_return, iter=>per_infected_sim)
        push!(to_return_abs, iter=>n_infected)

        push!(new_infected_direct_to_return, iter=>new_infected_direct)
        push!(new_infected_indirect_to_return, iter=>new_infected_indirect)
        push!(new_infected_both_to_return, iter=>new_infected_both)

        push!(susceptible_to_return, iter => susceptible)
    end
    

    Dict{Symbol, Any}(
        :infected_percentage => to_return,
        :n_infected => to_return_abs,
        :new_infected_direct => new_infected_direct_to_return,
        :new_infected_indirect => new_infected_indirect_to_return,
        :new_infected_both => new_infected_both_to_return,
        :susceptible_percentage => susceptible_to_return
    )
end
