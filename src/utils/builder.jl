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
        if get(loc2he, string(checkin.venueid), -1) == -1
            println(checkin.venueid)
            @assert(get(loc2he, string(checkin.venueid), -1) != -1)
        end

        # if a user visits the same place in the same timeframe
        # only the timestamp of his/her last checkin is stored
        setindex!(
            h,
            Dates.value(checkin.timestamp),  # checkin to store
            get(user2vertex, string(checkin.userid), -1), # node id
            get(loc2he, string(checkin.venueid), -1) # hyperedge id
        )
    end

    h
end