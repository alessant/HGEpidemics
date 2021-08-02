using HGEpidemics
using DataFrames
using Dates
using CSV
using JSON
using JSON3
using JSONTables
using Serialization

"""
   Squeeze check-ins
"""

############################
# Loading simulation params
############################
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
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf[1, :Δ])),
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )


smay = DateTime(2012,5,7)
Dates.monthname(smay)

function adjust_date(x)
    smay = DateTime(2012,5,7)
    emay = DateTime(2012,6,7)

    if x < smay
        return x + Dates.Month(1)
    elseif x > emay
        v = nothing

        if Dates.monthname(x) == "June"
            v = x - Dates.Month(1)
        elseif Dates.monthname(x) == "July"
            v = x - Dates.Month(2)
        elseif Dates.monthname(x) == "August"
            v = x - Dates.Month(3)
        else
            return x - Dates.Month(4)
        end
    else
        return x
    end
end

transformed_df = transform(df, :timestamp => x -> adjust_date.(x))

minimum(df[!, :timestamp])
minimum(transformed_df[!, :timestamp_function])

maximum(df[!, :timestamp])
maximum(transformed_df[!, :timestamp_function])

rename!(transformed_df, [:timestamp => :old_timestamp, :timestamp_function => :timestamp])
select!(transformed_df, Not(:old_timestamp))
CSV.write("data/foursquare/dataset_TSMC2014_TKY_JCS4_.txt", transformed_df)