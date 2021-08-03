using HGEpidemics
using DataFrames
using Dates
using CSV
using JSON
using JSON3
using JSONTables
using Serialization
using PyPlot

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
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf[1, :Δ])),
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )

april = (_start = Date(2012,4,7), _end = DateTime(2012,5,7))
may = (_start = Date(2012,5,7), _end = DateTime(2012,6,7))
june = (_start = Date(2012,6,7), _end = DateTime(2012,7,7))
july = (_start = Date(2012,7,7), _end = DateTime(2012,8,7))

Dates.monthname(april._start)
Dates.day(april._start)

months = [april, may, june, july]
month_names = ["April", "May", "June", "July"]

checkins_per_month = Dict{String, Dict{Int, Int}}()

for m in months
    _df = filter(
        r-> (r[:timestamp] >= m._start) &&
            (r[:timestamp] <= m._end),
            df
    )

    println(nrow(_df))

    month_data = get!(checkins_per_month, Dates.monthname(m._start), Dict{Int, Int}(collect(1:31) .=> 0))

    for checkin in eachrow(_df)
        # _monthname = Dates.monthname(checkin.timestamp)
        _day = Dates.day(checkin.timestamp)

        push!(
            month_data,
            _day => get!(month_data, _day, 0) + 1
        )
    end
end

checkins_per_month

dists_to_store = Vector{Int}[]

clf()

for month in month_names
    dist = [checkins_per_month[month][day] for day in 1:31]

    push!(dists_to_store, dist)

    plot(dist)
end

gcf()

sum(dists_to_store)
dists_to_store

open("$(output_path)/latex_n_checkins.csv", "w+") do f
    write(f, "day,all,april,may,june,july\n")

    for day in 1:31
        a = checkins_per_month["April"][day]
        m = checkins_per_month["May"][day]
        j = checkins_per_month["June"][day]
        jl = checkins_per_month["July"][day]

        all = a + m + j + jl

        to_write = "$day,$all,$a,$m,$j,$jl\n"
        
        write(f, to_write)
    end
end

