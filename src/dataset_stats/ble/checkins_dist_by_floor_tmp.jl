using Pkg
Pkg.activate(".")
using TVHEpidemicDynamics
using DataFrames
using Dates
using CSV
using JSON3
using PyPlot
using Statistics
using StatsBase
using PyPlot
using Serialization

"""
   Policy of social distancing
"""
project_path = dirname(pathof(TVHEpidemicDynamics))

output_path = joinpath(project_path, "experiments", "immunization", "ble", "results")


fdata_params =
    joinpath(project_path, "experiments", "spreading", "ble", "configs", "blebluetooth_dataset.json")


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

data_params.start_date

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
            Δ = convert(Dates.Millisecond, Dates.Hour(4)),
            δ = convert(Dates.Millisecond, Dates.Minute(15)),
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
serialize("loc_by_checkins.data", to_store)


#singoli intervalli
checkins = Array{Int, 1}()

checkins_I = Array{Int, 1}()
checkins_II = Array{Int, 1}()
checkins_III = Array{Int, 1}()

first_floor = ["00" "01" "02" "30" "25" "29" "08" "13" "28" "24" "23" "27" "26" "22"] #"31"
second_floor = ["14" "15" "16" "17" "18" "19" "20" "21"]
thrid_floor = ["03" "04" "05" "06" "07" "08" "09" "10" "11" "12"]

first_floor_ids = ["rpi-$(id)" for id in first_floor]
second_floor_ids = ["rpi-$(id)" for id in second_floor]
thrid_floor_ids = ["rpi-$(id)" for id in thrid_floor]


to_store_for_lockdown = filter(
    row -> (row.venueid in thrid_floor_ids),
    df
)

i = 0
for row in eachrow(df)
    if row.venueid in thrid_floor_ids
        global i += 1
    end
end
i


CSV.write("ble_closed_rooms.csv", to_store_for_lockdown, dateformat="Y-m-dTH:M:S")




to_store = [loc2he[id] for id in first_floor_ids]
serialize("ble_closed_rooms_I.data", to_store)


for i in 1:length(intervals)

    #key = Dates.hour(round(intervals[i].first, Dates.Hour))
    first, second, third = 0, 0, 0

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    push!(
        checkins,
        nrow(_df)
    )

    for c in eachrow(_df)
        if c.venueid in first_floor_ids
            first += 1
        elseif c.venueid in second_floor_ids
            second += 1
        else
            third += 1
        end
    end

    push!(checkins_I, first)
    push!(checkins_II, second)
    push!(checkins_III, third)
end




data = [checkins, checkins_I, checkins_II, checkins_III]


clf()
#figure(figsize=(7,4))
fig, ax1 = plt.subplots()

#ax1.set_ylim(bottom=0, top=5000)#, top=0.6)

for dist in data
    ax1.plot(dist)
end

ax1.legend(["all", "I", "II", "III"])

gcf()
savefig("checkins_ble.png")













#WEEK
first_floor = ["00" "01" "02" "30" "31" "25" "29" "08" "13" "28" "24" "23" "27" "26" "22"]
second_floor = ["14" "15" "16" "17" "18" "19" "20" "21"]
thrid_floor = ["03" "04" "05" "06" "07" "08" "09" "10" "11" "12"]

first_floor_ids = ["rpi-$(id)" for id in first_floor]
second_floor_ids = ["rpi-$(id)" for id in second_floor]
thrid_floor_ids = ["rpi-$(id)" for id in thrid_floor]

weekly_checkins = Dict{String, Array{Int, 1}}()

weekly_checkins_I = Dict{String, Array{Int, 1}}()
weekly_checkins_II = Dict{String, Array{Int, 1}}()
weekly_checkins_III = Dict{String, Array{Int, 1}}()


for i in 1:length(intervals)
    println(round(intervals[i].first, Dates.Hour))

    #key = Dates.hour(round(intervals[i].first, Dates.Hour))
    key = "$(Dates.dayname(intervals[i].first))_$(Dates.hour(round(intervals[i].first, Dates.Hour)))"
    first, second, third = 0, 0, 0

    # _df = filter(
    #     r-> (r[:timestamp] >= intervals[i].first) &&
    #         (r[:timestamp] <= intervals[i].second),
    #         df
    # )

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] < intervals[i].second),
            df
    )

    push!(
        get!(weekly_checkins, key, Array{Int, 1}()),
        nrow(_df)
    )

    for c in eachrow(_df)
        if c.venueid in first_floor_ids
            first += 1
        elseif c.venueid in second_floor_ids
            second += 1
        else
            third += 1
        end
    end

    push!(
        get!(weekly_checkins_I, key, Array{Int, 1}()),
        first
    )

    push!(
        get!(weekly_checkins_II, key, Array{Int, 1}()),
        second
    )

    push!(
        get!(weekly_checkins_III, key, Array{Int, 1}()),
        third
    )
end


weekly_checkins

mean(collect(values(weekly_checkins)))

mean.(values(weekly_checkins))

listed_data = [weekly_checkins, weekly_checkins_I, weekly_checkins_II, weekly_checkins_III]


days = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
abbr = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
hours = ["0","4","8", "12", "16","20"]

labels_dict = collect(zip(days, abbr))
labels_dict[1][2]

first = []
second = []
third = []

all = []
labels = []

for day in days
    for i in hours
        k = string(day, "_", i)
        push!(all, mean(weekly_checkins[k]))

        label = i == "0" ? string(day) : string(i)

        push!(labels, label)
    end    
end

map(elem -> (elem == "4" || elem == "8") ? string("0", elem) : elem , labels)
new_labels = map(elem -> (elem != "12" && !(elem in days)) ? "" : elem , labels)

new_labels = map(elem -> elem in days ? labels_dict[findall(day -> day == elem, days)[1]][2] : elem, new_labels)


findall(day -> day == "Monday", days)[1]

for day in days
    for i in hours
        k = string(day, "_", i)
        push!(first, mean(weekly_checkins_I[k]))
    end    
end

for day in days
    for i in hours
        k = string(day, "_", i)
        push!(second, mean(weekly_checkins_II[k]))
    end    
end

for day in days
    for i in hours
        k = string(day, "_", i)
        push!(third, mean(weekly_checkins_III[k]))
    end    
end



to_plot = [all, first, second, third]
linestyles = ["solid", "dashed", "dashdot", "dotted"]

clf()
figure(figsize=(7,4))
#fig, ax1 = plt.subplots()

ylim(bottom=0, top=3500)#, top=0.6)

for (index, dist) in enumerate(to_plot)
    plot(dist, linestyle=linestyles[index])
end

xticks(0:length(all)-1, labels=new_labels, rotation=0, fontweight="roman")
ylabel("Number of check-ins", fontweight="roman", fontsize="large", labelpad=10)

legend(["Entire building", "I-floor", "II-floor", "III-floor"])

#xticks(range(1, stop=35), step=0.2)
gcf()
tight_layout()
savefig("weekly_checkins_ble_tmis.png")




# DAY
first_floor = ["00" "01" "02" "30" "31" "25" "29" "08" "13" "28" "24" "23" "27" "26" "22"]
second_floor = ["14" "15" "16" "17" "18" "19" "20" "21"]
thrid_floor = ["03" "04" "05" "06" "07" "08" "09" "10" "11" "12"]

first_floor_ids = ["rpi-$(id)" for id in first_floor]
second_floor_ids = ["rpi-$(id)" for id in second_floor]
thrid_floor_ids = ["rpi-$(id)" for id in thrid_floor]


weekly_checkins = Dict{String, Array{Int, 1}}()

weekly_checkins_I = Dict{String, Array{Int, 1}}()
weekly_checkins_II = Dict{String, Array{Int, 1}}()
weekly_checkins_III = Dict{String, Array{Int, 1}}()


for i in 1:length(intervals)
    println(Dates.hour(round(intervals[i].first, Dates.Hour)))

    #key = Dates.hour(round(intervals[i].first, Dates.Hour))
    key = "$(Dates.hour(round(intervals[i].first, Dates.Hour)))"
    first, second, third = 0, 0, 0

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    push!(
        get!(weekly_checkins, key, Array{Int, 1}()),
        nrow(_df)
    )

    for c in eachrow(_df)
        if c.venueid in first_floor_ids
            first += 1
        elseif c.venueid in second_floor_ids
            second += 1
        else
            third += 1
        end
    end

    push!(
        get!(weekly_checkins_I, key, Array{Int, 1}()),
        first
    )

    push!(
        get!(weekly_checkins_II, key, Array{Int, 1}()),
        second
    )

    push!(
        get!(weekly_checkins_III, key, Array{Int, 1}()),
        third
    )
end


listed_data = [weekly_checkins, weekly_checkins_I, weekly_checkins_II, weekly_checkins_III]
hours = ["0", "4", "8", "12", "16", "20"]

first = []
second = []
third = []

all = []

to_plot = [all, first, second, third]
labels = []


for (index, value) in enumerate(listed_data)
    for h in hours
        push!(to_plot[index], mean(value[h]))

        label = (h == "0" || h == "4" || h == "8") ? string("0", h, "h00") : string(h, "h00") 
        push!(labels, label)
    end    
end

labels = ["00-04", "04-08", "08-12", "12-16", "16-20", "20-24"]
patterns = [ ".", "/" , "-", "x", "\\" ,  "|"  , "+" , "o", "O", "." ]

clf()
figure(figsize=(7,4))
#fig, ax1 = plt.subplots()

#ax1.set_ylim(bottom=0, top=3500)#, top=0.6)
x = 0:length(all)-1
width = 0.2
val = -1

for (index, dist) in enumerate(to_plot)
    global val 

    bar(x .+ width .* val, dist, width, alpha=0.8, hatch=patterns[index])
    val += 1
end

xticks(x .+ width .- .1, labels=labels, rotation=0)
ylabel(" ", fontweight="roman", fontsize="large", labelpad=10)
#xlim(collect(minimum(x)) .- width .* 2, collect(maximum(x)) .+ width*8 .- 0.4)

legend(["Entire building", "I-floor", "II-floor", "III-floor"])

#xticks(range(1, stop=35), step=0.2)
gcf()
tight_layout()
savefig("daily_checkins_ble.png")






# number of peole checking-in within one day
n_people = Array{Int, 1}()

for i in 1:length(intervals)

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    users = unique(df, Symbol(data_params.userid))

    push!(n_people, nrow(users))
end

most_crowded_day = findmax(n_people)
intervals[most_crowded_day[2]]

# we will only allow 31 people
# to access the building
df_cp = copy(df)
checkins_to_remove = Set{Int}()

for i in 1:length(intervals)
    added = 0

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    users = unique(df, Symbol(data_params.userid))
    chosen = sample(users[!, :userid], 31)

    for checkin in eachrow(_df)
        if !(checkin.userid in chosen)
            push!(checkins_to_remove, checkin.entry_id)
        end

    end
end

to_store = filter(
        row -> !(row.entry_id in checkins_to_remove),
        df
    )

CSV.write("ble.csv", to_store, dateformat="Y-m-dTH:M:S")

CSV.read("ble.csv", dateformat="Y-m-dTH:M:S")










#
# CLOSING ROOMS
#
order = ["00" "01" "02" "30" "31" "25" "29" "08" "13" "28" "24" "23" "27" "26" "22" "14" "15" "16" "17" "18" "19" "20" "21" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"]
_ids = ["rpi-$(id)" for id in order]

#rooms_to_close = _ids[1:2:length(_ids)]
rooms_to_close = ["30" "01" "29" "13" "24" "27" "22" "20" "18" "15" "5" "7" "11"]
_ids_rooms_to_close = ["rpi-$(id)" for id in rooms_to_close]

to_store_rooms = filter(
    row -> !(row.venueid in _ids_rooms_to_close),
    df
)

CSV.write("ble_closed_rooms.csv", to_store_rooms, dateformat="Y-m-dTH:M:S")


for row in eachrow(df)
    if row.venueid in rooms_to_close
end




###
# interval => people to verify whether they are infected or not
people_to_check = Dict{DateTime, Array{String, 1}}()

for i in 1:length(intervals)

    _df = filter(
        r-> (r[:timestamp] >= intervals[i].first) &&
            (r[:timestamp] <= intervals[i].second),
            df
    )

    users = unique(df, Symbol(data_params.userid))[!, Symbol(data_params.userid)]
    #chosen = sample(users[!, :userid], 10)
    users_string = [string(u) for u in users]

    push!(people_to_check, intervals[i].first => users_string)
end

using Serialization
serialize("people_to_check.data", people_to_check)

_df = filter(
    r-> (r[:timestamp] >= intervals[1].first) &&
        (r[:timestamp] <= intervals[1].second),
        df
)

user_groups = groupby(_df, Symbol(data_params.userid))

data = sort!(_df, [order(:userid), :timestamp])
data[data[:userid] .== 0, :]


using Serialization
deserialize("data/blebeacon/people_to_check.data")