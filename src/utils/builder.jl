"""
    generatehg(
        h::Hypergraph,
        df::DataFrame,
        mindate::DateTime,
        maxdate::DateTime,
        user2vertex::Dict{String, Int},
        loc2he::Dict{String, Int},
        usersepoc::Vector{Int},
        t::Int;
        all::Bool = true
    )

Populate a hypergraph `h` with checkin data from the dataframe `df`
happening within `mindate` and `maxdate`.

At the first step of the simulation, `h` is `nothing`. In this case,
it initializes a new hypergraph `h`, where its vertices correspond to users in `user2vertex`
and its hyperedges to locations in `loc2he`. For each user, the timestamp of the last checkin
he/she did in a given location is stored.

In the following simulation steps, `h` is modified according to the new checkin data.
Also in this case, the timestamp of the last checkin a user did in a given location is stored.

`usersepoc` containes information about whether a given user was
already present in the simulation. Based on this vector, the number of
newcoming users or users who just moved in another location is evaluated.

`t` is the current simulation step. In combination with `all`, it can be used
for debugging purposes.
"""
function generatehg!(
    h::Union{Hypergraph, Nothing},
    df::DataFrame,
    mindate::DateTime,
    maxdate::DateTime,
    user2vertex::Dict{String, Int},
    loc2he::Dict{String, Int},
    usersepoc::Vector{Int},
    t::Int;
    thr::Union{Int, Nothing} = nothing,
    all::Bool = true
)
    added, moved = 0, 0

    # select only the current timeframe
    currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

    # initialize hg
    if isnothing(h)
        # first step of the simulation
        # initialize the hypergraph from scratch
        h = inithg(currdf, user2vertex, loc2he; thr)
    else
        # clean the timestamp stored within the hypergraph
        # it stores only the last checkin for each user
        for v = 1:nhv(h)
            max, id = 0, 0
            for he in gethyperedges(h, v)
                if getindex(h, v, he.first) > max
                    max = getindex(h, v, he.first)
                    id = he.first
                end
                setindex!(h, nothing, v, he.first)
            end
            if id != 0
                setindex!(h, max, v, id)
            end
        end
    end

    # Repopulate the hypergraph with checkins
    # from the new timeframe
    for checkin in eachrow(currdf)
        he = get(loc2he, string(checkin.venueid), -1)
        users_per_location = length(getvertices(h, he))

        # we do not have to do anything here
        # we already have thr user in this location
        !isnothing(thr) && users_per_location == thr && continue
        
        # we already have more users than allowed
        # in this location
        # we need to remove the surplus
        if !isnothing(thr) && users_per_location > thr
            vs = [v.first for v in getvertices(h, he)]
            to_remove = fill(0, (users_per_location - thr))
            
            sample!(vs, to_remove; replace=false)

            for v in to_remove
                setindex!(h, nothing, v, he)
            end

            #println(users_per_location)
            # println(vs)
            # println(to_remove)
            # println(length(getvertices(h, he)))
        else 
            # we do not have any check on the thr or we still have room to insert agents in the locations

            # if the user was not in any other venue, i.e. an hyperedge,
            # then it is a new coming node
            usersepoc[get(user2vertex, string(checkin.userid), -1)] == 0 ?
                usersepoc[get(user2vertex, string(checkin.userid), -1)] = 1 : moved += 1

            # if a user visits the same place in the same timeframe
            # only the timestamp of his/her last checkin is stored
            if t == 1 || all #just for debugging purposes
                setindex!(
                    h,
                    Dates.value(checkin.timestamp), # checkin to store
                    get(user2vertex, string(checkin.userid), -1), # node id
                    get(loc2he, string(checkin.venueid), -1) # hyperedge id
                )
            end
        end
    end

    # number of newcoming users in the current timeframe
    added = nrow(currdf) - moved

    h, added, moved
end



"""
    inithg(
        df::DataFrame,
        user2vertex::Dict{String, Int},
        loc2he::Dict{String, Int}
    )

Initialize a hypergraph `h`, where its vertices correspond to users in `user2vertex`
and its hyperedges to locations in `loc2he`. For each user, the timestamp of the last checkin
he/she did in a given location is stored.

"""
function inithg(df::DataFrame, user2vertex::Dict{String, Int}, loc2he::Dict{String, Int}; thr::Union{Int, Nothing} = nothing)
    # Generate a hypergraph `h` with
    # `n` users and `m` locations
    # The ids of each user and each locations
    # are stored as metadata in `h`
    h = Hypergraph{Int, String, String}(length(keys(user2vertex)), length(keys(loc2he)))

    # Storing metadata
    for user in user2vertex
        set_vertex_meta!(h, user.first, user.second)
    end

    for place in loc2he
        set_hyperedge_meta!(h, place.first, place.second)
    end

    for checkin in eachrow(df)
        !isnothing(thr) && length(getvertices(h, get(loc2he, string(checkin.venueid), -1))) >= thr && continue

        # this should never happen
        # the df is cleaned from missing data
        if get(loc2he, checkin.venueid, -1) == -1
            println(checkin.venueid)
        end

        # if a user visits the same place in the same timeframe
        # only the timestamp of his/her last checkin is stored
        setindex!(
            h,
            Dates.value(checkin.timestamp),  # checkin to store
            get(user2vertex, string(checkin.userid), -1), # node id
            get(loc2he, checkin.venueid, -1) # hyperedge id
        )
    end

    h
end


# TODO adjust the following consider

function generate_hg_g_weighted!(
    df::DataFrame,
    mindate::Union{DateTime, Nothing},
    maxdate::DateTime,
    user2vertex::Dict{String, Int},
    loc2he::Dict{String, Int},
    δ::Dates.Millisecond
)
    currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

    h = Hypergraph{Int, String, String}(length(keys(user2vertex)), length(keys(loc2he)))
    hw = Hypergraph{Float64, String, String}(length(keys(user2vertex)), length(keys(loc2he)))
    g = SimpleGraph{Int}(length(keys(user2vertex)))
    gw = SimpleWeightedGraph(length(keys(user2vertex)))

    #map each node to a list of checkin in locid with time 
    # locid -> nodeid, time ; ....
    checkins = Dict{Int, Dict{Int, Int}}()
   
   
    for checkin in eachrow(currdf)
        # this should never happen
        # the df is cleaned from missing data
        if get(loc2he, checkin.venueid, -1) == -1
            println(checkin.venueid)
        end

        nodeid =  get(user2vertex, string(checkin.userid), -1)
        locid =  get(loc2he, checkin.venueid, -1)

        push!(
            get!(checkins, locid,  Dict{Int, Int}()),
            nodeid =>  Dates.value(checkin.timestamp)
        ) 
    end

    for locid in keys(checkins)
        visiting_nodes = get!(checkins, locid, nothing)
        isnothing(visiting_nodes) && continue

        for u in keys(visiting_nodes)
            ndirect = 0
            for v in keys(visiting_nodes)
                if u!=v && abs(visiting_nodes[u] -  visiting_nodes[v]) <= δ.value
                    ndirect += 1
                    add_edge!(g, u, v)
                    
                    w = gw.weights[u,v]
                    add_edge!(gw, u, v, w + 1.0)
                end
            end
            setindex!(
                h,
                1, #ndirect,  # number of direct contact of u in locid
                u, # node id
                locid # hyperedge id
            )

            setindex!(
                hw,
                ndirect/198, #ndirect,  # number of direct contact of u in locid
                u, # node id
                locid # hyperedge id
            )
        end
    end

    # for v in 1:nhv(h)
    #     for he in 1:nhe(h)
    #         if isnothing(h[v, he])
    #             h[v, he] = 0
    #         end
    #     end
    # end

    g, gw, h, hw
end


function evaluate_entropy(
    df::DataFrame,
    mindate::DateTime,
    maxdate::DateTime,
    user2vertex::Dict{String, Int},
    loc2he::Dict{String, Int}
)

    user2location = Dict{Int, Array{Int, 1}}()

    for user in user2vertex
        push!(
            user2location,
            user2vertex[user.first] => Array{Int, 1}()
        )
    end

    to_return = Dict{Int, Float64}()


    # for t in slots

    #     t == imm_start && continue

    #     mindate = get(intervals, t, 0).first
    #     maxdate = get(intervals, t, 0).second

    #     currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

    #     for checkin in eachrow(currdf)
    #         push!(
    #             get!(user2location, user2vertex[string(checkin.userid)], Array{Int, 1}()),
    #             loc2he[checkin.venueid]
    #         )
    #     end
    # end

    currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

    for checkin in eachrow(currdf)
        push!(
            get!(user2location, user2vertex[string(checkin.userid)], Array{Int, 1}()),
            loc2he[checkin.venueid]
        )
    end

    for u in keys(user2location)
        data = user2location[u]
        
        if length(data) == 0 || length(data) == 1 || length(unique(data)) == 1
            push!(
                to_return,
                u => 0 #typemin(Int)
            )
        else
            entropy = get_entropy(data)

            push!(
                to_return,
                u => entropy < eps() ? 0.0 : entropy
            )
        end
    end

    to_return
end



function evaluate_entropy_direct(
    df::DataFrame,
    mindate::DateTime,
    maxdate::DateTime,
    δ,
    user2vertex::Dict{String, Int},
    loc2he::Dict{String, Int}
)

    user2contact = Dict{Int, Array{Int, 1}}()
    to_return = Dict{Int, Float64}()


    #map each node to a list of checkin in locid with time 
    # locid -> nodeid, time ; ....
    checkins = Dict{Int, Dict{Int, Int}}()


    currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

    for checkin in eachrow(currdf)
        # this should never happen
        # the df is cleaned from missing data
        if get(loc2he, checkin.venueid, -1) == -1
            println(checkin.venueid)
        end

        nodeid =  get(user2vertex, string(checkin.userid), -1)
        locid =  get(loc2he, checkin.venueid, -1)

        push!(
            get!(checkins, locid,  Dict{Int, Int}()),
            nodeid =>  Dates.value(checkin.timestamp)
        ) 
    end


    for locid in keys(checkins)
        visiting_nodes = get!(checkins, locid, nothing)
        isnothing(visiting_nodes) && continue

        for u in keys(visiting_nodes)
            ndirect = 0
            for v in keys(visiting_nodes)
                if u!=v && abs(visiting_nodes[u] -  visiting_nodes[v]) <= δ.value
                    ndirect += 1
                    #add_edge!(g, u, v)

                    push!(  
                        get!(user2contact, v, Array{Int, 1}()),
                        u
                    )

                    push!(  
                        get!(user2contact, u, Array{Int, 1}()),
                        v
                    )
                end
            end
        end
    end


    # for t in slots

    #     t == imm_start && continue

    #     mindate = get(intervals, t, 0).first
    #     maxdate = get(intervals, t, 0).second

    #     currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

    #     g = build_graph(currdf, δ, user2vertex, loc2he)

    #     for v in vertices(g)
    #         for u in neighbors(g, v)
    #             push!(  
    #                 get!(user2contact, v, Array{Int, 1}()),
    #                 u
    #             )
    #         end
    #     end
    # end

    #currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)
    #g = build_graph(currdf, δ, user2vertex, loc2he)

    # for v in vertices(g)
    #     for u in neighbors(g, v)
    #         push!(  
    #             get!(user2contact, v, Array{Int, 1}()),
    #             u
    #         )
    #     end
    # end

    for u in keys(user2contact)
        data = user2contact[u]
        
        if length(data) == 1 || length(unique(data)) == 1
            push!(
                to_return,
                u => 0#typemin(Int)
            )
        else
            entropy = get_entropy(data)

            push!(
                to_return,
                u => entropy < eps() ? 0.0 : entropy
            )
        end
    end

    to_return
end



function evaluate_entropy_direct_mean(
    df::DataFrame,
    intervals,
    δ,
    user2vertex::Dict{String, Int},
    loc2he::Dict{String, Int},
    start_slot,
    end_slot
)

    entropies = Dict{Int, Array{Float64, 1}}()
    to_return = Dict{Int, Float64}()

    for t in start_slot:end_slot
        #
        user2contact = Dict{Int, Array{Int, 1}}()

        mindate = get(intervals, t, 0).first
        maxdate = get(intervals, t, 0).second

        currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

        g = build_graph(currdf, δ, user2vertex, loc2he)

        for v in vertices(g)
            for u in neighbors(g, v)
                push!(  
                    get!(user2contact, v, Array{Int, 1}()),
                    u
                )
            end
        end

        # evaluate evaluate_entropy
        for u in keys(user2contact)
            data = user2contact[u]
            
            if length(data) == 1
                push!(
                    get!(entropies, u, Array{Float64, 1}()),
                    0
                )
            else
                push!(
                    get!(entropies, u, Array{Float64, 1}()),
                    get_entropy(data)
                )
            end
        end
    end

    # for v in keys(entropies)
    #     λ = -1

    #     entropies_v = entropies[v]
    #     scaled_entropies = [entropy * ℯ ^ - λ * t for (t, entropy) in enumerate(entropies_v)]
    #     val = maximum(scaled_entropies)

    #     push!(to_return, v => val)
    # end

    for v in keys(entropies)
        push!(to_return, v => sum(entropies[v])/length(entropies[v]))
    end

    to_return
end




function evaluate_entropy_indirect_mean(
    df::DataFrame,
    intervals,
    user2vertex::Dict{String, Int},
    loc2he::Dict{String, Int},
    start_slot,
    end_slot
)

    entropies = Dict{Int, Array{Float64, 1}}()
    to_return = Dict{Int, Float64}()

    for t in range(start_slot, stop=end_slot)
        user2location = Dict{Int, Array{Int, 1}}()

        mindate = get(intervals, t, 0).first
        maxdate = get(intervals, t, 0).second

        currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

        for checkin in eachrow(currdf)
            push!(
                get!(user2location, user2vertex[string(checkin.userid)], Array{Int, 1}()),
                loc2he[checkin.venueid]
            )
        end

        #evaluate_entropy
        for u in keys(user2location)
            data = user2location[u]
            
            if length(data) == 1 || length(unique(data)) == 1
                push!(
                    get!(entropies, u, Array{Float64, 1}()),
                    0
                )
            else
                push!(
                    get!(entropies, u, Array{Float64, 1}()),
                    get_entropy(data)
                )
            end
        end
    end

    # for v in keys(entropies)
    #     λ = -1

    #     entropies_v = entropies[v]
    #     scaled_entropies = [entropy * ℯ ^ - λ * t for (t, entropy) in enumerate(entropies_v)]
    #     val = maximum(scaled_entropies)

    #     push!(to_return, v => val)
    # end
    
    for v in keys(entropies)        
        push!(to_return, v => sum(entropies[v]) / length(entropies[v]))
    end

    to_return
end




function evaluate_entropy_both(
    df::DataFrame,
    intervals,
    δ,
    user2vertex::Dict{String, Int},
    loc2he::Dict{String, Int},
)

    user2contact = Dict{Int, Array{Int, 1}}()
    to_return = Dict{Int, Float64}()

    for t in 1:length(intervals)

        t == 10 && continue

        mindate = get(intervals, t, 0).first
        maxdate = get(intervals, t, 0).second

        # indirect
        currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

        for checkin in eachrow(currdf)
            push!(
                get!(user2contact, user2vertex[string(checkin.userid)], Array{Int, 1}()),
                loc2he[checkin.venueid] + length(user2vertex)
            )
        end

        # direct
        g = build_graph(currdf, δ, user2vertex, loc2he)

        for v in vertices(g)
            for u in neighbors(g, v)
                push!(  
                    get!(user2contact, v, Array{Int, 1}()),
                    u
                )
            end
        end
    end

    for u in keys(user2contact)
        data = user2contact[u]
        
        if length(data) == 1
            push!(
                to_return,
                u => typemin(Int)
            )
        else
            push!(
                to_return,
                u => get_entropy(data)
            )
        end
    end

    to_return
end



function build_graph(df, δ, user2vertex, loc2he)
    g = SimpleGraph{Int}(length(keys(user2vertex)))

    #map each node to a list of checkin in locid with time 
    # locid -> nodeid, time ; ....
    checkins = Dict{Int, Dict{Int, Int}}()

    for checkin in eachrow(df)
        # this should never happen
        # the df is cleaned from missing data
        if get(loc2he, checkin.venueid, -1) == -1
            println(checkin.venueid)
        end

        nodeid =  get(user2vertex, string(checkin.userid), -1)
        locid =  get(loc2he, checkin.venueid, -1)

        push!(
            get!(checkins, locid,  Dict{Int, Int}()),
            nodeid =>  Dates.value(checkin.timestamp)
        ) 
    end


    for locid in keys(checkins)
        visiting_nodes = get!(checkins, locid, nothing)
        isnothing(visiting_nodes) && continue

        for u in keys(visiting_nodes)
            ndirect = 0
            for v in keys(visiting_nodes)
                if u!=v && abs(visiting_nodes[u] -  visiting_nodes[v]) <= δ.value
                    ndirect += 1
                    add_edge!(g, u, v)
                end
            end
        end
    end

    g
end

