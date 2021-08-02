using Distances
using HGEpidemics
using NSGAII
using PyPlot
using Serialization

path = "src/experiments/NSGA/FQ/results/both/lock-crowded/15_iter/"
name = "res-100-100"

data = deserialize("$(path)$(name).data")

##############################
# CONFIGS
##############################

# MIN INFECTED
min_infected_data = sort(data, by = x -> x.y[1])[1]
min_infected_conf = min_infected_data.pheno
min_infected_point = round.(min_infected_data.y, digits=2)
HGEpidemics.parse_conf(min_infected_conf)


# min_damage_data = sort(data, by = x -> x.y[2])[1]
min_damage_data = sort(data, by = x -> (x.y[2] + x.y[3]))[1]
min_damage_conf = min_damage_data.pheno
min_damage_point = round.(min_damage_data.y, digits=2)
HGEpidemics.parse_conf(min_damage_conf)

# min_distance_data = sort(data, by = x -> euclidean([0, 0], [x.y[1], x.y[2]]))[1]
min_distance_data = sort(data, by = x -> euclidean([0, 0, 0], [x.y[1], x.y[2], x.y[3]]))[1]
min_distance_conf = min_distance_data.pheno
min_distance_point = round.(min_distance_data.y, digits=2)
HGEpidemics.parse_conf(min_distance_conf)


##############################
# STORE DATA FOR SPIDER CHART
##############################
name = "agents_spider"
output_path = "$(path)paper/$(name).data"

HGEpidemics.save_spider_data(output_path, min_infected_conf, min_damage_conf, min_distance_conf)


##############################
# PLOT
##############################
pop = filter(indiv -> indiv.rank <= 1, data) #keeps only the non-dominated solutions

min_infected_point_color = "blue"
min_damage_point_color = "red"
min_distance_point_color = "green"

colors = Array{String, 1}()

for x in pop
    point = round.(x.y, digits=2)

    if point == min_infected_point
        push!(colors, min_infected_point_color)
    elseif point == min_damage_point
        push!(colors, min_damage_point_color)
    elseif point == min_distance_point
        push!(colors, min_distance_point_color)
    else
        push!(colors, "grey")
    end
end
colors

clf()
# scatter(map(x -> x.y[1], data), map(x -> x.y[2], data), c = colors) 
scatter(map(x -> x.y[1], pop), map(x -> x.y[2], pop), c = colors) 

# scatter3D(map(x -> x.y[1], data), map(x -> x.y[2], data), map(x -> x.y[3], data), c = colors)

xlabel("Infected nodes in %", fontweight="semibold", labelpad=10, fontsize="x-large")
ylabel("Damage over \n agents in %", fontweight="semibold", fontsize="x-large", labelpad=10)
# zlabel("Damage over \n locations in %", fontweight="semibold", fontsize="x-large", labelpad=10)
        
tight_layout()
gcf()
savefig("$(path)pareto-both.png")


##############################
# STORE DATA FOR SCATTER PLOT
##############################
name = "agents_scatter"
output_path = "$(path)/paper/$(name).data"

pop = filter(indiv -> indiv.rank <= 1, data) #keeps only the non-dominated solutions

min_infected_point_class = "min-infected"
min_damage_point_class = "min-damage"
min_distance_point_class = "min-distance"

classes = Array{String, 1}()

for x in pop
    point = round.(x.y, digits=2)

    if point == min_infected_point
        push!(classes, min_infected_point_class)
    elseif point == min_damage_point
        push!(classes, min_damage_point_class)
    elseif point == min_distance_point
        push!(classes, min_distance_point_class)
    else
        push!(classes, "other")
    end
end
classes

HGEpidemics.save_scatter_data(output_path, pop, classes)


