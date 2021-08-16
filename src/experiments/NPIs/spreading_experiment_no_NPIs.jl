using Pkg
Pkg.activate(".")

using HGEpidemics

using CSV
using DataFrames
using Dates
using JSON
using JSON3
using JSONTables
using PyPlot
using Statistics

"""

    Experiments on a spreading process
    via Time-Varying Hypergraphs
    using the Susceptible-Infected-Susceptible epidemic model. 

"""

############################
# Loading simulation params
############################
path = ARGS[1] 
input_data = JSON.parse((open(path, "r")))

output_path = input_data["output_path"]
fdata_params = input_data["data_params"]
fparams = input_data["sim_params"]

#project_path = dirname(pathof(HGEpidemics))

jtable = jsontable(read(open(fparams, "r")))
paramsdf = DataFrame(jtable)

# just a trick to group together
# all experiments to show in the same plot
test_data = Dict{String, Array{Any, 1}}()
for params in eachrow(paramsdf)
    push!(
        get!(test_data, params[:exp_id], Array{Any, 1}()),
        params
    )
end


#########################
# Generating model data
########################
data_params = JSON3.read(read(open(fdata_params, "r")))
header = [Symbol(col) for col in data_params.header]

# The choice of the interval within which
# either an indirect (Δ) or direct (δ) contact
# may occur influences the data the
# simulation is run on.
# For this reason, it is necessary to store
# diffent information according to the
# values of both Δ and δ.
intervals = unique(paramsdf, [:Δ, :δ])[!, [:Δ, :δ]]
intervals_data = Dict{String, Dict{Symbol, Any}}()

# df_ok = CSV.read(
#     data_params.dataset,
#     DataFrame;
#     datarow = 2,
#     header = header,
#     dateformat = data_params.dateformat,
#     limit = 10
# )

# println(data_params.dataset)

for i in eachrow(intervals)
    df, _intervals, user2vertex, loc2he =
        generate_model_data(
            data_params.dataset,
            header,
            Symbol(data_params.userid),
            Symbol(data_params.venueid),
            Symbol(data_params.UTCtime),
            data_params.dateformat;
            datarow = data_params.datarow,
            Δ = convert(Dates.Millisecond, Dates.Hour(i.Δ)),
            δ = convert(Dates.Millisecond, Dates.Minute(i.δ)),
            maxdate = Dates.DateTime(data_params.end_date),
            mindate = Dates.DateTime(data_params.start_date)
        )

    push!(
        get!(intervals_data, "$(i.Δ)$(i.δ)", Dict{Symbol, Any}()),
        :df => df,
        :intervals => _intervals,
        :user2vertex => user2vertex,
        :loc2he => loc2he,
    )
end



#########################
# Initialization of infected nodes
########################

# For the reproducibility of the experiments,
# the infected nodes at start as to be the same
per_infected = unique(paramsdf, [:infected_percentage])[!, [:infected_percentage]]
per_infected_data = Dict{Float64, Array{Int, 1}}()

users = keys(intervals_data[collect(keys(intervals_data))[1]][:user2vertex])

for p in eachrow(per_infected)
    vstatus = fill(1, length(users))
    vrand = rand(Float64, (1, length(users)))

    for i=1:length(users)
        if p.infected_percentage  <= vrand[i]
            vstatus[i] = 0
        end
    end

    vstatus = fill(0, length(users))
    vstatus[1] = 1

    push!(
        per_infected_data,
        p.infected_percentage => vec(vstatus)
    )
end


#########################
# Simulation
########################
simulation_data = Dict{String, Array{Pair{String, NamedTuple}, 1}}()

for testtype in keys(test_data)
    for test in get(test_data, testtype, nothing)

        println("----------------EXP CONFIG-------------------------")
        for property in propertynames(test)
            print("$(property) = $(test[property])  |   ")
        end
        println("\n---------------------------------------------------")

        runningparams = get(intervals_data, "$(test[:Δ])$(test[:δ])", Dict{Symbol, Any}())

        res_path =
            joinpath(output_path, "csv", "$(test[:exp_id])_$(test[:exp])_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).csv")

        results =
            simulate(
                SIS_NPIs(),
                get!(runningparams, :df, nothing),
                get!(runningparams, :intervals, nothing),
                get!(runningparams, :user2vertex, nothing),
                get!(runningparams, :loc2he, nothing),
                convert(Dates.Millisecond, Dates.Minute(test[:δ]));
                Δ = test[:Δ],
                vstatus = per_infected_data[test[:infected_percentage]],
                per_infected = test[:infected_percentage],
                c = test[:c],
                βd = test[:βd], 
                βᵢ = test[:βᵢ],
                βₑ = test[:βₑ],
                γₑ = test[:γₑ],
                γₐ = test[:γₐ],
                niter = 10,
                output_path = res_path,
                store_me = false
            )

        # get the average over all iterations
        infected_distribution = mean(collect(values(results[:infected_percentage])))
        susceptible_distribution = mean(collect(values(results[:susceptible_percentage])))
        
        agents_met = mean(collect(values(results[:agents_met])))
        location_visited = mean(collect(values(results[:location_visited])))
    
        push!(
            get!(simulation_data, testtype, Array{Dict{String, NamedTuple}, 1}()),
            test[:label] => (
                infected_distribution = infected_distribution, 
                infected_percentage_per_run = results[:infected_percentage],
                susceptible_distribution = susceptible_distribution,
                susceptible_percentage_per_run = results[:susceptible_percentage],
                agents_met = agents_met,
                location_visited = location_visited,
                Δ = test[:Δ], 
                δ = test[:δ]
            )
        )
    end
end


#########################
# Plotting infected ditribution
########################
if input_data["plot_infected"]
    plot_infected_distribution(default_plot(),simulation_data; output_path=output_path)
end

if haskey(input_data, "store_infected_distribution") && input_data["store_infected_distribution"]
    store_infected_distribution_data(simulation_data; output_path=output_path)
end


#########################
# store sim data
########################
if input_data["store_simulation_data"]
    using Serialization

    my_title = "simulation_data_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).data"
    println("Saving data in ... $(output_path)/$(my_title)")
    serialize("$output_path/$my_title", simulation_data)
end

#########################
# SIS data
########################
if haskey(input_data, "plot_SIS") && input_data["plot_SIS"]
    for k in keys(simulation_data)
        infected_per_run = simulation_data[k][1].second.infected_percentage_per_run # direct/indirect
        susceptible_per_run = simulation_data[k][1].second.susceptible_percentage_per_run
        
        plot_SIS_distribution(infected_per_run, susceptible_per_run; output_path=output_path)
    end
end

if haskey(input_data, "store_SIS_data") && input_data["store_SIS_data"]
    for k in keys(simulation_data)
        infected_per_run = simulation_data[k][1].second.infected_percentage_per_run # direct/indirect
        susceptible_per_run = simulation_data[k][1].second.susceptible_percentage_per_run
        
        store_SIS_distribution_data(infected_per_run, susceptible_per_run; output_path=output_path)
    end
end


#########################
# Storing damage
########################
if haskey(input_data, "store_damage") && input_data["store_damage"]
    timeframe = collect(keys(intervals_data))[1]
    user2vertex = intervals_data[timeframe][:user2vertex]

    test_name = collect(keys(simulation_data))[1]
    sim_data = simulation_data[test_name][1].second.agents_met
    sim_data_loc = simulation_data[test_name][1].second.location_visited

    damage_output_path_agents = input_data["damage_output_path_agents"]
    damage_output_path_locs = input_data["damage_output_path_locs"] 

    store_damage(sim_data, user2vertex, damage_output_path_agents)
    store_damage(sim_data_loc, user2vertex, damage_output_path_locs)
end


