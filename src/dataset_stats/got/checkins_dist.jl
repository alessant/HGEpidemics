using Pkg
Pkg.activate(".")

using HGEpidemics
using DataFrames
using Dates
using CSV
using JSON
using JSON3
using JSONTables
using PyPlot
using Serialization
using StatsBase
using Statistics

############################
# Loading simulation params
############################
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
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf[1, :Δ])),
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )



    
# number of locations
unique(df, :venueid)

characters = unique(df[!, :userid])

# building data for (#direct, #locations) scatter PyPlot
user2direct = Dict{String, Int}()
user2locations = Dict{String, Int}()

# building data for unique(#direct, #locations) scatter PyPlot
user2direct_unique = Dict{String, Int}()
user2locations_unique = Dict{String, Int}()

user2contacts = Dict{String, Set{String}}()

for c in characters
    push!(
        user2contacts,
        c => Set{String}()
    )
end

# evaluate the number of different locations per character
groups_by_char = groupby(df, :userid)

for char in groups_by_char
    nlocs = nrow(char) 
    nlocs_unique = nrow(unique(char, :venueid))
    user = char[1, :userid]
    
    push!(user2locations, user => nlocs)
    push!(user2locations_unique, user => nlocs_unique)
end
user2locations
user2locations_unique

# check
sort(collect((user2locations_unique)), by = x -> x[2], rev=true)

# evaluate the number of direct contacts
groups_by_scene = groupby(df, :timestamp)

for scene in groups_by_scene
    chars = scene[!, :userid]
    contacts = length(chars) - 1

    for c in chars
        push!(
            user2direct,
            c => get!(user2direct, c, 0) + contacts
        )

        for c1 in chars
            length(chars) == 1 && continue

            if c != c1
                push!(
                    get!(user2contacts, c, Set{String}()),
                    c1
                )
            end
        end
    end    
end
user2direct
user2contacts

for char in keys(user2contacts)
    push!(
        user2direct_unique,
        char => length(user2contacts[char])
    )
end
user2direct_unique

# check
sort(collect((user2direct_unique)), by = x -> x[2], rev=true)

##################################
# BASIC PLOT
##################################
characters = keys(user2direct)
x = collect(values(user2direct))
y = collect(values(user2locations))

cor(x, y)

clf()
scatter(x, y)
gcf()

# outliers
sort(collect(user2direct), by = x -> x[2], rev=true)

clf()
hist(x)
gcf()

##################################
# BASIC PLOT UNIQUE
##################################
characters = keys(user2direct)
x = collect(values(user2direct_unique))
y = collect(values(user2locations_unique))

cor(x, y)

clf()
scatter(x, y)
gcf()


open("$(output_path)/got_data.csv", "w+") do fp
    println(fp, "character,num_direct,num_loc")

    for c in characters
        to_write = string(c, ",", user2direct[c], ",", user2locations[c])
        println(fp, to_write)
    end
end

open("$(output_path)/got_data_unique.csv", "w+") do fp
    println(fp, "character,num_direct,num_loc")

    for c in characters
        to_write = string(c, ",", user2direct_unique[c], ",", user2locations_unique[c])
        println(fp, to_write)
    end
end