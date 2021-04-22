using HGEpidemics

using CSV
using Dates
using DataFrames
using JSON3
using JSONTables
using PyPlot


"""
 This script analyzes, over the most crowded week of the data set,
 the following elements:

 - The daily distribution of the differences (in seconds) between
 two consecutive check-ins within the same location;

 - The daily check-in distribution within the same place
 considering only direct contact;

 - The daily check-in distribution within the same place
 considering only indirect contact. This corresponds to evaluate
 the distribution of how many different locations each user visits per day.
"""

############################
# Loading data set params
############################
project_path = dirname(pathof(HGEpidemics))

# FOURSQUARE
output_path = joinpath(project_path, "dataset_stats", "foursquare", "plots")
fparams = joinpath(project_path, "dataset_stats", "foursquare", "configs", "foursquare.json")


############################
# Reading data
############################
data_params = JSON3.read(read(open(fparams, "r")))
header = [Symbol(col) for col in data_params.header]

df = CSV.read(
    data_params.dataset,
    DataFrame;
    copycols = true,
    header = header,
    dateformat = data_params.dateformat
    )

df = dropmissing(df)

# getting the first checkin within the data
mindate = minimum(df[!, :UTCtime])
# we evaluate the data starting from the midnight of the first monday
first_mon = floor(Dates.tonext(d -> Dates.dayofweek(d) == Dates.Monday, mindate), Day)


############################
# Evaluating the weekly check-in distribution
############################

# Find the weeks to analyze
Δₕ = 24*7
Δₘ = convert(Dates.Millisecond, Dates.Hour(Δₕ))

# interval_id -> (start_date, end_date)
intervals = find_intervals(first_mon, maximum(df[!, :UTCtime]), Δₘ)

# Count checkin data according to the discrete time intervals
# evaluated at the previous step
#
# interval_id -> # of check-ins
checkins_per_interval = evaluate_checkin_density(intervals, df)

# Sorting weeks according to their number of checkins
sorted_checkins_per_interval =
    sort(collect(checkins_per_interval), by = x -> x[2], rev = true)

# Working with the week containing
# the majority of checkins
most_crowded_week = intervals[sorted_checkins_per_interval[1].first]



"""
 Analyzing:
 - check-ins distribution;
 - direct contacts distribution;
 - indirect contact distribution.

 Here, we will use a discretization interval of δₕ = 4 hours.
"""
δₕ = 4
δₘ = convert(Dates.Millisecond, Dates.Hour(δₕ))

# interval_id -> (start_date, end_date)
intervals_δ = find_intervals(most_crowded_week.first, most_crowded_week.second, δₘ)

# CHECKINS DISTRIBUTION
# Evaluate the distribution of the differences (in seconds) between
# two consecutive check-ins within the same location.
distr = evaluate_checkins_distribution(intervals_δ, df)
plt_name = "checkins_distribution_most_crowded_week.png"

# DIRECT CONTACTS DISTRIBUTION
# Evaluatin direct contacts distribution within each interval.
distr = evaluate_direct_contacts_distribution(intervals_δ, df, convert(Dates.Millisecond, Dates.Hour(1)))
plt_name = "direct_contacts_distribution_most_crowded_week.png"

# INDIRECT CONTACTS DISTRIBUTION
# Evaluating indirect contact distribution (i.e. location distribution)
distr = evaluate_location_distribution(intervals_δ, df)
plt_name = "location_distribution_most_crowded_week.png"



"""
    Plot checkin distribution within the most crowded week:
    - Foursquare: 7th-13th May
"""
clf()

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)

ylabel("Time in Seconds", fontweight="semibold")
xlabel("Time Intervals", fontweight="semibold")

for t=1:length(keys(distr))
    values = get(distr, t, Array{Any, 1}())
    ax.boxplot(values, widths=0.5, positions=[t])
end

ax.set_xlim(0, 43)

# foursquare dates
xlabels = ["7 May", "04", "08", "12", "16", "20", "8 May", "04", "08", "12", "16", "20", "9 May", "04", "08", "12", "16", "20", "10 May", "04", "08", "12", "16", "20", "11 May", "04", "08", "12", "16", "20", "12 May", "04", "08", "12", "16", "20", "13 May", "04", "08", "12", "16", "20"]

xticks(1:42, xlabels, rotation=80)

gcf()

plt.tight_layout(.5)

savefig(joinpath(output_path, plt_name))