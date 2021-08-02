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
            datarow = data_params.datarow,
            Δ = convert(Dates.Millisecond, Dates.Hour(paramsdf[1, :Δ])),
            δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ])),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )


#
#
#
δ = convert(Dates.Millisecond, Dates.Minute(paramsdf[1, :δ]))
data_direct = Dict{Int, Array{Int, 1}}()
data_indirect = Dict{Int, Array{Int, 1}}()

for t=1:length(intervals)
    mindate = intervals[t].first
    maxdate = intervals[t].second

    currdf = filter(r->((r.timestamp >= mindate) && (r.timestamp < maxdate)), df)

    groups_by_user = groupby(currdf, :userid)

    # INDIRECT
    # evaluate the number of different locations per user
    if length(groups_by_user) == 0
        push!(get!(data_indirect, t, Array{Int, 1}()), 0)
    else
        for user_df in groups_by_user
            nlocs_unique = nrow(unique(user_df, :venueid))
            user = user_df[1, :userid]
            
            push!(
                get!(data_indirect, t, Array{Int, 1}()),
                nlocs_unique
            )
        end
    end

    # DIRECT
    # evaluate the number of direct contacts 
    contacts = Dict{Union{Int,String}, Set{Union{Int,String}}}()
    groups_by_loc = groupby(currdf, :venueid)

    #! assumption: already ordered
    for loc_df in groups_by_loc
        for checkin₁ in eachrow(loc_df) #! 1:n
            for checkin₂ in eachrow(loc_df) #! j = i+1:n
                "entry_id" in names(currdf) && checkin₁.entry_id == checkin₂.entry_id && continue
                checkin₁.userid == checkin₂.userid && continue

                if abs(checkin₁.timestamp - checkin₂.timestamp) <= δ
                    user₁ = checkin₁.userid
                    user₂ = checkin₂.userid

                    push!(get!(contacts, user₁, Set{Int}()), user₂)
                    push!(get!(contacts, user₂, Set{Int}()), user₁)
                end
                #! else continue
            end
        end
    end

    println(t)

    push!(
        data_direct,
        t => collect(length.(values(contacts)))
    )  
    
end
data_direct
data_indirect



to_plot_direct = [Array{Int, 1}(), Array{Float64, 1}(), Array{Float64, 1}()]
to_plot_indirect = [Array{Int, 1}(), Array{Float64, 1}(), Array{Float64, 1}()]


# DATA TO PLOT
for t=1:length(intervals)

    if length(data_direct[t]) == 0
        #median
        push!(to_plot_direct[1], 0)
        # 25%
        push!(to_plot_direct[2], 0)
        # 75%
        push!(to_plot_direct[3], 0)
    else
        q_direct = quantile(data_direct[t], [0.0, 0.25, 0.5, 0.75, 1.0])

        push!(to_plot_direct[1], q_direct[3])
        push!(to_plot_direct[2], q_direct[2])
        push!(to_plot_direct[3], q_direct[4])
    end

    if length(data_indirect[t]) == 0
        #median
        push!(to_plot_indirect[1], 0)
        # 25%
        push!(to_plot_indirect[2], 0)
        # 75%
        push!(to_plot_indirect[3], 0)
    else
        q_indirect = quantile(data_indirect[t], [0.0, 0.25, 0.5, 0.75, 1.0])

        #median
        push!(to_plot_indirect[1], q_indirect[3])
        # 25%
        push!(to_plot_indirect[2], q_indirect[2])
        # 75%
        push!(to_plot_indirect[3], q_indirect[4])
    end
end

for elem in to_plot_indirect
    println(length(elem))
end

# PLOT
clf()
figure(figsize=(17,7))

#direct
plot(to_plot_direct[1], marker="x", markevery=5, markersize=7)
fill_between(
    range(0, stop=length(intervals)-1), 
    to_plot_direct[2], 
    to_plot_direct[3],
    alpha=0.3
    )

#INDIRECT
plot(to_plot_indirect[1], marker="<", markevery=5, markersize=7)
fill_between(
    range(0, stop=length(intervals)-1), 
    to_plot_indirect[2], 
    to_plot_indirect[3],
    alpha=0.3
    )

tick_params(labelsize="large")
ylabel("Number of contacts", fontweight="roman", fontsize="x-large", labelpad=10)
xlabel(L"\Delta t", fontweight="roman", fontsize="x-large", labelpad=10)

legend(["Direct", "Indirect"])

title("GoT complete")

tight_layout(.5)
gcf()


open("$(output_path)/latex-contact-dist.csv", "w+") do f
    write(f, "x,direct,direct-25,direct-75,indirect,indirect-25,indirect-75\n")

    for i in 1:length(to_plot_direct[1])
        direct = to_plot_direct[1][i]
        direct_25 = to_plot_direct[2][i]
        direct_75 = to_plot_direct[3][i]

        indirect = to_plot_indirect[1][i]
        indirect_25 = to_plot_indirect[2][i] 
        indirect_75 = to_plot_indirect[3][i]

        write(
            f,
            "$i,$direct,$direct_25,$direct_75,$indirect,$indirect_25,$indirect_75\n"
        )
    end
end


savefig("$(output_path)/contact_dist.png") 