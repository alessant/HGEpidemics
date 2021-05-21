module HGEpidemics

using CSV
using DataFrames
using Dates
using Random
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


include("dataset_stats/utils.jl")

include("epidemics/sim_types.jl")
include("epidemics/SIS.jl")
include("epidemics/SIS_NPIs.jl")

include("epidemics/selection_strategies.jl")

include("epidemics/sim_utils.jl")
include("epidemics/logging_utils.jl")

include("utils/loader.jl")
include("utils/builder.jl")

end # module
