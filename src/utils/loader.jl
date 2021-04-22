"""
    generate_model_data()

Load a dataframe from a file `fname` with columns `columns`, containing checkin-data.

Note: the name of the columns have to be Symbol.

Return:
- a `dataframe` containing all the checkins that will be used for the simulation;
- a `dict` of intervals within which an indirect contact may occur:
    - at each interval `i` corresponds a pair of dates delimiting that interval
    - `i` => (date₁, date₂)
- a `dict` maping each user to the corresponding node id in the hypergraph
- a `dict` maping each location to the corresponding hyperedge id in the hypergraph

### Arguments
- `userid` : column corresponding to user ids (and to vertices);
- `venueid` : column corresponding to venue ids (and to hyperedges);
- `timestamp` : column corresponding to the checkin timestamps (and to weights);
- `dateformat`: the format of the timestamp of the checkin;
- `Δ`: time within which an indirect contact may happen;
- `δ`: time within which a direct contact may happen;
- `limit`: number of rows to read. Each row represents a checkin;
- `mindate`: if set, a checkin is considered only if it happend after `mindate`;
- `maxdate`: if set, a checkin is considered only if it happend before `maxdate`.
"""
function generate_model_data(
    fname::AbstractString,
    columns::AbstractVector{Symbol},
    userid::Symbol,
    venueid::Symbol,
    timestamp::Symbol,
    dateformat::AbstractString;
    Δ::Millisecond = Dates.Millisecond(14400000),
    δ::Millisecond = Dates.Millisecond(1000),
    limit::Integer = typemax(Int32),
    mindate::Union{DateTime,Nothing} = nothing,
    maxdate::Union{DateTime,Nothing} = nothing
)

    df = CSV.read(
        fname,
        DataFrame;
        header = columns,
        dateformat = dateformat,
        #limit = limit
    )

    # remove all rows without a timestamp
    df = dropmissing(df, timestamp)
    df = dropmissing(df)

    # filter out all the rows not in the
    # timeframe to consider
    df = filter(
        r-> (isnothing(mindate) || r[timestamp] >= mindate) &&
            (isnothing(maxdate) || r[timestamp] <= maxdate),
            df
    )

    sort!(df, [timestamp])

    # vertices
    users = unique(df, userid)
    numusers = nrow(users)

    # hyperedges - upper buond
    places = unique(df, venueid)
    numplaces = nrow(places)

    # time slots
    mintime = minimum(df[!, timestamp])
    maxtime = maximum(df[!, timestamp])

    nintervals = ceil(Int, (maxtime - mintime) / Δ)
    intervals = Dict{Int, Pair{DateTime, DateTime}}()

    # build timeframes
    for i = 1:nintervals
        offset =
            mintime + Δ + δ > maxtime ? maxtime + Dates.Millisecond(1) : mintime + Δ + δ

        get!(
            intervals,
            i,
            Pair{DateTime,DateTime}(
                convert(Dates.DateTime, (mintime - δ)), convert(Dates.DateTime, (offset))
            )
        )

        mintime = offset - δ
    end

    # mapping user -> node
    user2vertex = Dict{String, Int}()
    index = 1

    # add vertices
    for user in eachrow(users)
        uid = string(user[userid])
        if !haskey(user2vertex, uid)
            push!(user2vertex, uid => index)
            index += 1
        end
    end

    # mapping place -> he
    loc2he = Dict{String, Int}()
    index = 1

    # add he
    for place in eachrow(places)
        placeid = string(place[venueid])
        if !haskey(loc2he, placeid)
            push!(loc2he, placeid => index)
            index += 1
        end
    end

    rename!(df, [userid => :userid, venueid => :venueid, timestamp => :timestamp])

    df, intervals, user2vertex, loc2he
end
