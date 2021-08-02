using Pkg
Pkg.activate(".")

using HGEpidemics

using Distributions 
using JSON
using NSGAII
using Plots
using Random
using Serialization

Random.seed!(10);

# GLOBAL PARAMS - SIMULATION DATA
data_params = ARGS[1] #"src/experiments/NSGA/BLE/configs/ble_params.json"
println("Working with $(data_params)")

println("Init simulation params...")
simulation_data = init_simulation(data_params)

data_params_dict = JSON.parse((open(data_params, "r")))

# caching results
individuals = Dict{Any, NamedTuple}()


####################################
# UTILS
####################################
function discretize_individual(individual)
    new_individual = [round(npi, digits=2) for npi in individual]
    new_individual
end

function create_population(n_pop)
    init_pop = [rand(Uniform(0.0,1.0), 7) for _ in 1:n_pop]
    return init_pop
end


####################################
# NSGA - simulation
####################################
function simulate(individual)
    global individuals
    global simulation_data
    global data_params_dict

    d_individual = discretize_individual(individual)

    if !haskey(individuals, hash(d_individual))
        simulation_res = run_simulation(individual, simulation_data)

        # saving fitness
        infected = simulation_res.infected_distribution[length(simulation_res.infected_distribution)]
        agent_damage = check_damage(simulation_res.agents_met, data_params_dict["agent_damage_path"])
        location_damage = check_damage(simulation_res.location_visited, data_params_dict["location_damage_path"])

        to_return = (infected = infected, agent_damage = agent_damage, location_damage = location_damage)

        # caching res
        push!(
            individuals,
            hash(d_individual) => to_return
        )

        return to_return
    else
        individuals[hash(d_individual)]
    end
end

####################################
# NSGA - objective functions
####################################
function obj_function_agents(x) 
    simulate(x).infected, simulate(x).agent_damage
end

function obj_function_locations(x) 
    simulate(x).infected, simulate(x).location_damage
end

function obj_function_both(x) 
    simulate(x).infected, simulate(x).agent_damage, simulate(x).location_damage
end


####################################
# plot
####################################
function plot_pop(pop, path)
    pop = filter(indiv -> indiv.rank <= 1, pop) #keeps only the non-dominated solutions
    scatter(map(x -> x.y[1], pop), map(x -> x.y[2], pop)) #, markersize = 1) #|> display
    savefig(path)
    #sleep(0.1)
end


####################################
# NSGA PARAMS
####################################
npop = 100 #200
ngen = 100 #100

#Encodes two variables 0 <= x_i <= 1, with a precision of 1E-4
const bc = BinaryCoding(4, [0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1]) 

# obj function
if data_params_dict["damage"] == "agents"
    obj_function = obj_function_agents
elseif data_params_dict["damage"] == "locations"
    obj_function = obj_function_locations
else
    obj_function = obj_function_both
end


initial_pop = create_population(npop)
#
# 1. popsize
# 2. ngen
# 3. evaluation function
# 4. initialization function
# seed = [[1.,-1.],[2.5,0.5],[0.5,0.25]]

println("Starting optimization...")
flush(stdout)

res = nsga(npop, ngen, obj_function, bc, seed = initial_pop) #, fplot = plot_pop)

path = data_params_dict["output_path"]
serialize(joinpath(path, "res-$npop-$ngen.data"), res)
plot_pop(res, joinpath(path, "res-$npop-$ngen.png"))
