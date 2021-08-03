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

"""
   Plot checkins distribution
"""
############################
# Loading simulation params
############################
path = "src/dataset_stats/ble/configs/ble_params.json"
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

    key = "$(Dates.dayname(intervals[i].first))_$(Dates.hour(round(intervals[i].first, Dates.Hour)))"
    first, second, third = 0, 0, 0

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

# findall(day -> day == "Monday", days)[1]

for day in days
    for i in hours
        k = string(day, "_", i)
        push!(first, mean(weekly_checkins_I[k]))
        push!(second, mean(weekly_checkins_II[k]))
        push!(third, mean(weekly_checkins_III[k]))
    end    
end

to_plot = [all, first, second, third]
linestyles = ["solid", "dashed", "dashdot", "dotted"]

clf()
figure(figsize=(7,4))

ylim(bottom=0, top=3500)

for (index, dist) in enumerate(to_plot)
    plot(dist, linestyle=linestyles[index])
end

xticks(0:length(all)-1, labels=new_labels, rotation=0, fontweight="roman")
ylabel("Number of check-ins", fontweight="roman", fontsize="large", labelpad=10)

legend(["Entire building", "I-floor", "II-floor", "III-floor"])

gcf()
tight_layout()


open("$(output_path)/latex-ble-weekly-checkins.csv", "w+") do f
    write(f, "all,first,second,third,x\n")

    for i in 1:length(all)
        to_write = "$(all[i]),$(first[i]),$(second[i]),$(third[i]),$i\n"
        write(f, to_write)
    end
end

open("$(output_path)/latex-ble-weekly-checkins-xtickslabels.csv", "w+") do f
    
    for i in 1:length(new_labels)
        write(f, "$(new_labels[i]),")
    end
end



# DAY
first_floor = ["00" "01" "02" "30" "31" "25" "29" "08" "13" "28" "24" "23" "27" "26" "22"]
second_floor = ["14" "15" "16" "17" "18" "19" "20" "21"]
thrid_floor = ["03" "04" "05" "06" "07" "08" "09" "10" "11" "12"]

first_floor_ids = ["rpi-$(id)" for id in first_floor]
second_floor_ids = ["rpi-$(id)" for id in second_floor]
third_floor_ids = ["rpi-$(id)" for id in thrid_floor]

weekly_checkins = Dict{String, Array{Int, 1}}()

weekly_checkins_I = Dict{String, Array{Int, 1}}()
weekly_checkins_II = Dict{String, Array{Int, 1}}()
weekly_checkins_III = Dict{String, Array{Int, 1}}()

for i in 1:length(intervals)
    println(Dates.hour(round(intervals[i].first, Dates.Hour)))

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

legend(["Entire building", "I-floor", "II-floor", "III-floor"])

gcf()
tight_layout()

open("$(output_path)/latex-ble-daily-checkins.csv", "w+") do f
    write(f, "x,all,first,second,third\n")

    for i in 1:length(all)
        to_write = "$(labels[i]),$(all[i]),$(first[i]),$(second[i]),$(third[i])\n"
        write(f, to_write)
    end
end


open("$(output_path)/latex-ble-daily-checkins-xtickslabels.csv", "w+") do f
    
    for i in 1:length(labels)
        write(f, "$(labels[i]),")
    end
end
