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

path = "src/dataset_stats/got/configs/got_params.json"
input_data = JSON.parse((open(path, "r")))

output_path = input_data["output_path"]
fdata_params = input_data["data_params"]
fparams = input_data["sim_params"]

df = CSV.read("$(output_path)/got_data.csv", DataFrame)
df = CSV.read("$(output_path)/got_data_unique.csv", DataFrame)

df[!, :index] .= range(1, stop=nrow(df))

#attach color
df[!, :color] .= "#0080ff"
df[df[!, :character] .== "Jon_Snow", :color] .= Ref("#e6ac00")
df[df[!, :character] .== "Daenerys_Targaryen", :color] .= Ref("#e6ac00")
df[df[!, :character] .== "Tyrion_Lannister", :color] .= Ref("#e6ac00")

clf()

x = df[!, :num_direct]
y = df[!, :num_loc]
colors = df[!, :color]

# left, width = 0.1, 0.65
# bottom, height = 0.1, 0.65
# spacing = 0.005

#rect_scatter = [left, bottom, width, height]
# rect_histx = [left, bottom + height + spacing, width, 0.2]
# rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(figsize=(7, 7))

#ax_scatter = plt.axes(rect_scatter)
tick_params(direction="in", top=true, right=true, labelsize="large")

xlabel(L"d_w(\mathcal{G})", fontweight="roman", fontsize="xx-large", labelpad=10)
ylabel(L"d_w(\mathcal{H})", fontweight="roman", fontsize="xx-large", labelpad=10)

# ax_histx = plt.axes(rect_histx)
# ax_histx.tick_params(direction="in", labelbottom=false)
# ax_histy = plt.axes(rect_histy)
# ax_histy.tick_params(direction="in", labelleft=false)

# the scatter plot:
scatter(x, y, marker="x", c=colors)

# ax_histx.hist(x, bins=50, color="#0080ff")
# ax_histy.hist(y, orientation="horizontal", bins=50, color="#0080ff")

for row in eachrow(filter(row -> row.character in ["Jon_Snow", "Daenerys_Targaryen", "Tyrion_Lannister"], df))
    name = split(row.character, "_")
    name = string(name[1], " ", name[2])

    offset_x = 1.2 
    offset_y = -6

    # offset_x = 1.2 
    # offset_y = -0.2 

    if name == "Jon Snow"
        offset_x = -510

        # offset_x = -30
        # offset_y = -0.3
    elseif name == "Daenerys Targaryen"
        name = "Daenerys\nTargaryen"
        offset_x = -250
        offset_y = -50

        # offset_x = -60 #-900
        # offset_y = -0.3
    elseif name == "Tyrion Lannister"
        offset_x = -830

        # name = "  Tyrion\nLannister"
        # offset_x = -13 #-900
        # offset_y = -2.5
    end

    annotate(name, (row.num_direct + offset_x, row.num_loc + offset_y), fontsize="x-large")
end

tight_layout(.5)
gcf()

#########################
# Save plot
########################
savefig("$(output_path)/degree_dist.png")
