using HGEpidemics
using DataFrames
using Dates
using CSV
using JSON
using JSON3
using JSONTables
using Serialization
using StatsBase

"""
   Pre-compute the most crowded locations
"""

############################
# Loading simulation params
############################
# path = "src/dataset_stats/ble/configs/ble_params.json"
path = "src/dataset_stats/foursquare/configs/fq_params.json"
input_data = JSON.parse((open(path, "r")))

output_path = input_data["output_path"]
fdata_params = input_data["data_params"]
fparams = input_data["sim_params"]

jtable = jsontable(read(open(fparams, "r")))
paramsdf = DataFrame(jtable)
data_params = JSON3.read(read(open(fdata_params, "r")))

# The choice of the interval within which
# either an indirect (Δ) or direct (δ) contact
# may occur influences the data the
# simulation is run on.
# For this reason, it is necessary to store
# diffent information according to the
# values of both Δ and δ.
intervals_data = Dict{String, Dict{Symbol, Any}}()

header = [Symbol(col) for col in data_params.header]

# evaluating new dataset
# removing people
df, intervals, user2vertex, loc2he =
        generate_model_data(
            data_params.dataset,
            header,
            Symbol(data_params.userid),
            Symbol(data_params.venueid),
            Symbol(data_params.UTCtime),
            data_params.dateformat;
            datarow = data_params.datarow,
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf[1, :Δ])),
            δ = convert(Dates.Millisecond, Dates.Minute(0)), #δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )

# evaluating new dataset
# removing people


function filter_checkins(df)
    agents = unique(df[:, :userid])
    n_agents = length(agents)

    #n_agents <= 3 && return df

    chosen = sample(agents, ceil(Int, n_agents/2))

    to_return = filter(
        row -> row.userid in chosen,
        df
    )

    to_return
end



new_df = DataFrame(
    userid = Int64[], 
    venueid = String[], 
    catid = String[],
    catname = String[], 
    lat = Float64[],
    lng = Float64[],
    timezoneoffset = Int64[],
    timestamp = DateTime[]
)

df_rows = 0
n_check = 0

for i in 1:length(intervals)
    global new_df, n_check, df_rows

    # filter df within Δ
    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    df_rows += nrow(_df)
    sum = 0

    # we are interested in reducing the access
    # to a place within Δ
    gdf = groupby(_df, :venueid)
    #loc_users_df = combine(gdf, :userid => length∘unique => :n_agents)

    for group in gdf
        filtered_df = filter_checkins(group)
        #println("$(nrow(_df)) -- $(nrow(group)) --- $(nrow(filtered_df))")

        #isnothing(filtered_df) && continue

        #println(filtered_df)
        # println(names(filtered_df))
        #append!(new_df, filtered_df; cols = :orderequal)
        append!(new_df, filtered_df; cols = :orderequal)

        sum += (nrow(filtered_df))
        n_check += (nrow(filtered_df))
    end

    println(nrow(_df), " - ", sum)
end
new_df
n_check
df_rows

CSV.write("data/foursquare/dataset_TSMC2014_NYC_social_distancing_JSC4.txt", new_df)

gdf = groupby(df, :venueid)
gdf[1]
combine(gdf, :userid => length∘unique => :n_agents)

for group in gdf
    #filtered_df = filter_checkins(group)
    println(nrow(group))# -- $(nrow(filtered_df))")
end

names(df_ok)