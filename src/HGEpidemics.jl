module HGEpidemics

using CSV
using DataFrames
using Dates
using JSON
using JSON3
using JSONTables
using PyPlot
using Random
using Revise
using Serialization
using SimpleHypergraphs
using StatsBase


export find_intervals, evaluate_checkin_density, evaluate_checkins_distribution
export evaluate_direct_contacts_distribution, evaluate_location_distribution

export Abstract_Simulation_Model
export SIS, SIS_NPIs
export simulate

export uniform
export custom

export f, infected, infected_ppm
export init_log

export generate_model_data
export inithg, generatehg!
export plot_infected_distribution, plot_SIS_distribution
export store_infected_distribution_data, store_SIS_distribution_data

export store_damage, check_damage
export init_simulation, run_simulation

include("dataset_stats/utils.jl")

include("epidemics/sim_types.jl")
include("epidemics/SIS.jl")
include("epidemics/SIS_NPIs.jl")

include("epidemics/selection_strategies.jl")

include("epidemics/sim_utils.jl")
include("epidemics/logging_utils.jl")

include("utils/loader.jl")
include("utils/builder.jl")
include("utils/plotter.jl")
include("utils/latex.jl")

include("experiments/NSGA/damage.jl")
include("experiments/NSGA/sim_utils.jl")

include("experiments/NSGA/latex_utils.jl")

end # module
