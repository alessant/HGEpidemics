using LaTeXStrings
using PyPlot
using Serialization
using Statistics
using TVHEpidemicDynamics

############à
# 5.3 time varying
###########
sim = deserialize("src/experiments/immunization/got/results/aamas21/data/got/exp53.data")
_keys = ["global_53_$(p)" for p in range(10, stop=100, step=10)]

n_examples = length(sim[_keys[1]][1].second.infected_distribution) 
n_examples_imm = n_examples - 70

# all techniques
#approaches = ["g_degree", "g_pr", "Direct_entropy", "hg_degree", "hg_pr", "Indirect_entropy"]
approaches = ["g_degree", "Direct_entropy", "hg_degree", "Indirect_entropy"]

# unmitigated
unmitigated_dist = sim[_keys[1]][1].second.infected_raw
unmitigated_mean = mean([v[n_examples] for v in values(unmitigated_dist)])
unmitigated_data = repeat([unmitigated_mean], length(_keys))

# uniform line
uniform_dist = sim[_keys[1]][2].second.infected_raw
val_mean = mean([v[n_examples] for v in values(uniform_dist)])
uniform_data = repeat([val_mean], length(_keys))


# other distributions
to_plot = Dict{String, Array{Float64, 1}}()
to_plot_std = Dict{String, Array{Float64, 1}}()

for index in 3:length(sim[_keys[1]])
    for exp in _keys
        dist = sim[exp][index].second.infected_raw
        data = [v[n_examples] for v in values(dist)]
    
        m = mean(data)
        s = std(data)
    
        push!(
            get!(to_plot, sim[exp][index].first, Array{Float64, 1}()),
            m
        )

        push!(
            get!(to_plot_std, sim[exp][index].first, Array{Float64, 1}()),
            s
        )
    end
end
to_plot
to_plot_std  



# PLOT0
labels = [L"No \quad intervention", L"Uniform", L"d(\mathcal{G})", L"\mathbb{H}[\mathcal{K}_\mathcal{D}]", L"d(\mathcal{H})", L"\mathbb{H}[\mathcal{K}_\mathcal{E}]"]
markers = ["o", "d", "v", "s"]

clf()
figure(figsize=(7,4))

plot(unmitigated_data, marker="x")
fill_between(range(0, stop=length(unmitigated_data)-1), unmitigated_data.-0, unmitigated_data.+0 ,alpha=0.3)

plot(uniform_data, marker="+", color="#763568")
fill_between(range(0, stop=length(uniform_data)-1), uniform_data.-0, uniform_data.+0 ,alpha=0.3)

for (index, app) in enumerate(approaches)
    data = to_plot[app]
    stds = to_plot_std[app]

    plot(data, marker=markers[index])
    fill_between(range(0, stop=length(data)-1), data-stds, data+stds ,alpha=0.3) #, facecolor=clrs[i]
    #errorbar(range(0, stop=(length(data)-1)), data, yerr=stds)
end

ylim(top=1)

legend(labels, ncol=3, loc="best", bbox_to_anchor=(0.5, 0., 0.5, 0.9))

tick_params(labelsize="large")
xticks(range(0, stop=length(unmitigated_data)-1), range(.1, stop=1, step=.1), fontweight="roman", fontsize="large")

ylabel("Infected (%)", fontweight="roman", fontsize="x-large", labelpad=10)
xlabel("Agent knowledge (%)", fontweight="roman", fontsize="x-large", labelpad=10)

tight_layout(.5)
gcf()
close()

#########################
# Save plot
########################
project_path = dirname(pathof(TVHEpidemicDynamics))
output_path = joinpath(project_path, "experiments", "immunization", "got", "results", "aamas21", "plot")

savefig("$(output_path)/exp53.png")