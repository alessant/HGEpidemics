module HGEpidemics

using CSV
using DataFrames
using Dates
using SimpleHypergraphs


export find_intervals, evaluate_checkin_density, evaluate_checkins_distribution
export evaluate_direct_contacts_distribution, evaluate_location_distribution

export Abstract_Simulation_Model
export SIS
export simulate

export f, infected
export init_log

export generate_model_data
export inithg, generatehg!


include("dataset_stats/utils.jl")

include("epidemics/sim_types.jl")
include("epidemics/SIS.jl")

include("epidemics/sim_utils.jl")
include("epidemics/logging_utils.jl")

include("utils/loader.jl")
include("utils/builder.jl")

end # module
