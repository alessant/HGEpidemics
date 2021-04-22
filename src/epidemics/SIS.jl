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
            γₐ::Float64 = 0.1
            niter::Int = 1,
            output_path::Union{AbstractString, Nothing} = nothing,
            print_me::Bool = true,
            store_me::Bool = true,
            kwargs...
    )

 Simulate a Susceptible-Infected-Susceptible diffusion model exploiting a
 Time-Varying Hypergraph. An immunization strategy may be applied either
 on vertices (i.e. people) and hyperedges (i.e. locations).

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
 - `niter`, number of iteration the simulation is repeated;
 - `output_path`, path where logs are stored;
 - `kwargs`, other optional params.

"""
function simulate(
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
        βd::Union{Int,Float64} = 0.2,
        βₑ::Union{Int,Float64} = 0.06,
        βᵢ::Union{Int,Float64} = 0.1,
        γₑ::Union{Int,Float64} = 0.06,
        γₐ::Union{Int,Float64} = 0.1,
        niter::Int = 10,
        output_path::Union{AbstractString, Nothing} = nothing,
        print_me::Bool = true,
        store_me::Bool = true,
        kwargs...
)

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
            "moved_agents,perc_infected_agents,perc_infected_locations\n"
            )

        init_log(output_path, header)
    end


    #########################
    # Init simulation params
    ########################

    # iter -> percentage of infected agents per simulation step
    to_return = Dict{Int, Array{Float64, 1}}()

    new_infected_direct_to_return = Dict{Int, Array{Int, 1}}()
    new_infected_indirect_to_return = Dict{Int, Array{Int, 1}}()
    new_infected_both_to_return = Dict{Int, Array{Int, 1}}()

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

        new_infected_direct = Array{Int, 1}()
        new_infected_indirect = Array{Int, 1}()
        new_infected_both = Array{Int, 1}()

        # store which agents are present in the given timeframe
        agentsepoc = zeros(Int, length(user2vertex))

        # Initialize the status vector of the nodes
        _vstatus = copy(vstatus)
        # Initially, all location are clean
        hestatus = zeros(Int, length(loc2he))

        # Storing the new status of each vertex and hyperedge
        vnextstatus = copy(vstatus)
        henextstatus = copy(hestatus)

        push!(per_infected_sim, sum(_vstatus) / length(_vstatus))


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
            # DIFFUSION ALGORITHM
            ########################

            #
            # PHASE 1 - Agent-to-Environment
            #
            for he=1:nhe(h)

                # If the location has at least two agents
                # and it is not infected, it may become contamined.
                if length(getvertices(h, he)) > 1 && hestatus[he] == 0
                    i = infected(h, he, _vstatus)
                    if rand() <  1 - ℯ ^ - (βₑ * f(i, c))
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
                                    if _vstatus[v] == 0
                                        i += _vstatus[u.first]
                                    end
                                    avg_direct_contacts += 1
                                end
                            end
                        end
                    end
                    # a agent becomes infected according to
                    # the following probability
                    if _vstatus[v] == 0 && rand() < 1 - ℯ ^ - (βd * i)
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

                    # if the agent is healthy 
                    # it may become infected
                    if _vstatus[v] == 0
                        i = 0

                        for he in gethyperedges(h, v)
                            if length(getvertices(h, he.first)) > 1
                                i += hestatus[he.first]
                            end
                        end

                        if rand() < 1 - ℯ ^ -(βᵢ * i) #1 - ℯ^-(βᵢ * f(i, c))
                            vnextstatus[v] = 1
                            new_infected_ind += 1
                        end

                    elseif rand() < 1 - ℯ ^ - γₐ
                        # the agent spontaneously returns healthy
                        vnextstatus[v] = 0
                    end
                end
            end

            _vstatus = copy(vnextstatus)
            hestatus = copy(henextstatus)

            push!(per_infected_sim, sum(_vstatus) / length(_vstatus))

            push!(new_infected_direct, new_infected_dir)
            push!(new_infected_indirect, new_infected_ind)
            push!(new_infected_both, (new_infected_dir+new_infected_ind))


            to_write = string(
                "$(t),$(Δ),$(δ),$(c),$(per_infected),$(βd),$(βₑ),$(βᵢ),$(γₑ),$(γₐ),",
                "$(avg_he_size),$(avg_degree),$(avg_direct_contacts),$(added),$(moved),",
                "$(sum(_vstatus)/length(_vstatus)),$(sum(hestatus)/length(hestatus))\n"
                )

            store_me && log(output_path, to_write)
        end

        push!(to_return, iter => per_infected_sim)

        push!(new_infected_direct_to_return, iter => new_infected_direct)
        push!(new_infected_indirect_to_return, iter => new_infected_indirect)
        push!(new_infected_both_to_return, iter => new_infected_both)
    end
    

    to_return, new_infected_direct_to_return, new_infected_indirect_to_return, new_infected_both_to_return
    Dict{Symbol, Any}(
        :infected_percentage => to_return,
        :new_infected_direct => new_infected_direct_to_return,
        :new_infected_indirect => new_infected_indirect_to_return,
        :new_infected_both => new_infected_both_to_return
    )
end


