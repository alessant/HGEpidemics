using LaTeXStrings
using PyPlot
using Serialization
using Statistics
using TVHEpidemicDynamics

sim = deserialize("src/experiments/immunization/got/results/aamas21/data/got/exp51.data")
sim_original = deserialize("src/experiments/immunization/got/results/aamas21/data/got/baseline_original.data")
key = "baseline"

n_examples = length(sim[key][1].second.infected_distribution) 

# [0-70]
origin_data = sim_original[key][1].second.infected_distribution

linestyles = ["solid", "solid", "solid", "dashed", "dashdot", "dotted"]
markers = ["", "<", "o", "", "", ""]

linestyle = 1
marker = 1

clf()
figure(figsize=(7,4))

#labels = ["Direct&Indirect", "Direct", "Indirect", L"\frac{Direct&Indirect}{2}", L"\frac{Direct&Indirect}{4}", "Direct&Indirect t=0"]
labels = [L"\mathcal{D}\mathcal{E}", L"\mathcal{D}", L"\mathcal{E}", L"\frac{\mathcal{D}\mathcal{E}}{2}", L"\frac{\mathcal{D}\mathcal{E}}{4}", L"\mathcal{D}\mathcal{E} \quad t=0"]
labels = [L"\mathcal{D}\mathcal{E}", L"\mathcal{D}", L"\mathcal{E}", L"{\mathcal{D}\mathcal{E}}/{2}", L"{\mathcal{D}\mathcal{E}}/{4}", L"\mathcal{D}\mathcal{E} \quad t=0"]

for (index, exp) in enumerate(sim[key])
    global linestyle, marker
    data = exp.second.infected_distribution

    plot(data, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=5)

    linestyle = (linestyle + 1) % (length(linestyles)+1)
    marker = (marker + 1) % (length(markers)+1)

    if linestyle == 0
        linestyle = 1
    end
    if marker == 0
        marker = 1
    end
end

plot(origin_data)

ylim(top=1)

tick_params(labelsize="large")

ylabel("Infected (%)", fontweight="roman", fontsize="x-large", labelpad=10)
xlabel("Time intervals", fontweight="roman", fontsize="x-large", labelpad=10)

legend(labels, ncol=1)#, fontsize="large")#, loc="best", bbox_to_anchor=(0.5, 0., 0.5, 0.9))

tight_layout(.5)
gcf()


#########################
# Save plot
########################
project_path = dirname(pathof(TVHEpidemicDynamics))
output_path = joinpath(project_path, "experiments", "immunization", "got", "results", "aamas21", "plot")

savefig("$(output_path)/exp51_new.png")



















