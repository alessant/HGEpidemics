module HGEpidemics

using CSV
using DataFrames
using Dates
using DataStructures
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
include("dataset_stats/utils.jl")


export Abstract_Simulation_Model
export SIS, SIS_NPIs, SIS_NPIs_vax
export simulate
include("epidemics/sim_types.jl")
include("epidemics/SIS.jl")
include("epidemics/SIS_NPIs.jl")
include("epidemics/SIS_NPIs_vax.jl")


export store_damage, check_damage
include("epidemics/damage.jl")


export uniform
export custom
include("epidemics/selection_strategies.jl")

export degrees
export centrality
export random_walk
export acquaintance
export hg_pagerank
export whg_pagerank
export g_pagerank
export wg_pagerank
export hg_entropy
export g_entropy
export hg_entropy_mean
export g_degrees
export g_entropy_mean
export both_entropy
export both_entropy_norm
export both_entropy_mean
export acquaintance_new
export most_loved
include("epidemics/immunization_strategies.jl")


export f, infected, infected_ppm
export init_log
include("epidemics/sim_utils.jl")
include("epidemics/logging_utils.jl")


export generate_model_data
export inithg, generatehg!
include("utils/loader.jl")
include("utils/builder.jl")


export initialize_params
include("utils/params_handler.jl")


export Abstract_Plot_Type
export default_plot, infected_not_isolated, infected_not_quarantined
export isolation, quarantine
export plot_infected_distribution, plot_status_distribution, plot_SIS_distribution
include("utils/plotter.jl")


export store_infected_distribution_data, store_SIS_distribution_data
include("utils/latex.jl")


export init_simulation, run_simulation
include("experiments/NSGA/sim_utils.jl")
include("experiments/NSGA/latex_utils.jl")

end # module
