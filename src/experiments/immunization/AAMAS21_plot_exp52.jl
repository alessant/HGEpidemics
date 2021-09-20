using LaTeXStrings
using PyPlot
using Serialization
using Statistics
using TVHEpidemicDynamics

############à
# 5.2 GRAFH vs HYPERGRAPH
###########
sim_52 = deserialize("src/experiments/immunization/got/results/aamas21/data/got/exp52_local.data")
key = "local_52"

n_examples = length(sim_52[key][1].second.infected_distribution) 
n_examples_imm = n_examples - 70

# all techniques
approaches = ["g_degree", "g_pr", "Direct_entropy", "hg_degree", "hg_pr", "Indirect_entropy"]
approaches = ["Acq_uniform", "Acq_degree_g", "Acq_direct_entropy", "Acq_degree_hg", "Acq_indirect_entropy"] #local

# unmitigated
unmitigated_dist = sim_52[key][1].second.infected_raw
unmitigated_mean = mean([v[n_examples] for v in values(unmitigated_dist)])
unmitigated_data = repeat([unmitigated_mean], length(approaches))

# uniform line
uniform_dist = sim_52[key][2].second.infected_raw
val_mean = mean([v[n_examples] for v in values(uniform_dist)])
uniform_data = repeat([val_mean], length(approaches))

# other distributions
to_plot = Dict{String, Array{Float64, 1}}()

for index in 3:length(sim_52[key])
    dist = sim_52[key][index].second.infected_raw
    data = [v[n_examples] for v in values(dist)]

    push!(
        to_plot,
        sim_52[key][index].first => data
    )
end
to_plot


# PLOT

data = [to_plot[app] for app in approaches]
labels = [L"d(\mathcal{G})", L"pr(\mathcal{G})", L"$\mathbb{H}[\mathcal{K}_\mathcal{D}]$", L"d(\mathcal{H})", L"pr(\mathcal{H})", L"$\mathbb{H}[\mathcal{K}_\mathcal{E}]$"]
labels = [L"acq. uniform", L"acq. d(\mathcal{G})", L"$acq. \mathbb{H}[\mathcal{K}_\mathcal{D}]$", L"acq. d(\mathcal{H})", L"$acq. \mathbb{H}[\mathcal{K}_\mathcal{E}]$"] #local

clf()

figure(figsize=(7,4))

plot(unmitigated_data, marker="x")
plot(uniform_data, marker="+", color="#763568")
boxplot(data, positions = range(0, stop=length(data)-1))

ylim(top=1) #bottom=0.0, 

tick_params(labelsize="large")
xticks(range(0, stop=length(data)-1), labels, fontweight="roman", fontsize="x-large")

ylabel("Infected (%)", fontweight="roman", fontsize="x-large", labelpad=10)

legend(["No intervention", "Uniform"], fontsize="large", loc="best", bbox_to_anchor=(0.5, 0., 0.5, 0.4))

tight_layout(.5)
gcf()
close()


#########################
# Save plot
########################
project_path = dirname(pathof(TVHEpidemicDynamics))
output_path = joinpath(project_path, "experiments", "immunization", "got", "results", "aamas21", "plot")

savefig("$(output_path)/exp52_local.png")