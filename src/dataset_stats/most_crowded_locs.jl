using Pkg
Pkg.activate(".")
using HGEpidemics
using DataFrames
using Dates
using CSV
using JSON
using JSON3
using JSONTables
using Serialization

"""
   Pre-compute the most crowded locations
"""

############################
# Loading simulation params
############################
# path = "src/dataset_stats/ble/configs/ble_params.json"
# path = "src/dataset_stats/foursquare/configs/fq_params.json"
path = "src/dataset_stats/got/configs/got_params.json"
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
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )

unique(df[!, :userid])
unique(df[!, :venueid])

# CHECK-INS DISTRIBUTION WITHIN ROOMS
checkins_loc = Dict{String, Int}()

for (key, subdf) in pairs(groupby(df, :venueid))
    println("Number of data points for $(key.venueid): $(nrow(subdf))")
    push!(checkins_loc, key.venueid => nrow(subdf))
end

sorted = sort(collect(checkins_loc), by = x -> x[2], rev = true)
to_store = [loc2he[id.first] for id in sorted]
# serialize("$(output_path)/loc_by_checkins.data", to_store)

loc_to_check = Set([sorted[elem].first for elem in 1:166])

_df = unique(df, :venueid)
categories = unique(filter(row -> row.venueid in loc_to_check, _df)[!, :catname])

categories