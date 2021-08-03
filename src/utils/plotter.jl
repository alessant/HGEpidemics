"""
    Plotting types
"""
abstract type Abstract_Plot_Type end

struct default_plot <: Abstract_Plot_Type end
struct infected_not_isolated <: Abstract_Plot_Type end
struct infected_not_quarantined <: Abstract_Plot_Type end

struct isolation <: Abstract_Plot_Type end
struct quarantine <: Abstract_Plot_Type end

"""
"""
function plot_infected_distribution(type::default_plot, simulation_data; output_path="")
    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]
    
    for test_type in keys(simulation_data)
        linestyle = 1
        marker = 1
        labels = Array{String, 1}()
        mytitle = "infected_dist_$(test_type)_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"
    
        clf()
        figure(figsize=(7,4))
    
        for exp in get!(simulation_data, test_type, Array{Float64, 1}())
            ylim(bottom=0.0, top=1.0)
            plot(exp.second.infected_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)
    
            xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
            ylabel("Δ = $(exp.second.Δ) hours \n Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)
            title("δ = $(exp.second.δ) minutes", pad=10, fontweight="semibold", fontsize="x-large")
    
            tick_params(labelsize="large")
    
            push!(labels, exp.first)
    
            linestyle = (linestyle + 1) % (length(linestyles)+1)
            marker = (marker + 1) % (length(markers)+1)
    
            if linestyle == 0
                linestyle = 1
            end
            if marker == 0
                marker = 1
            end
        end

        legend(labels, fontsize="large", ncol=2)
        plt.tight_layout()

        savefig("$(output_path)/$(mytitle)")
        println("Infected distribution -> saving figure in ... $(output_path)/$(mytitle)\n")
    end    
end


"""
"""
function plot_infected_distribution(type::infected_not_isolated, simulation_data; output_path="")
    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]
    
    for test_type in keys(simulation_data)
        linestyle = 1
        marker = 1
        labels = Array{String, 1}()
        mytitle = "infected_not_isolated_dist_$(test_type)_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"
    
        clf()
        figure(figsize=(7,4))
    
        for exp in get!(simulation_data, test_type, Array{Float64, 1}())
            ylim(bottom=0.0, top=1.0)
            plot(exp.second.infected_not_isolated_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)
    
            xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
            ylabel("Δ = $(exp.second.Δ) hours \n Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)
            title("δ = $(exp.second.δ) minutes", pad=10, fontweight="semibold", fontsize="x-large")
    
            tick_params(labelsize="large")
    
            push!(labels, exp.first)
    
            linestyle = (linestyle + 1) % (length(linestyles)+1)
            marker = (marker + 1) % (length(markers)+1)
    
            if linestyle == 0
                linestyle = 1
            end
            if marker == 0
                marker = 1
            end
        end
        legend(labels, fontsize="large", ncol=2)
        plt.tight_layout()
        savefig("$(output_path)/$(mytitle)")
        println("Infected not isolated distribution -> saving figure in ... $(output_path)/$(mytitle)\n")
    end    
end


"""
"""
function plot_infected_distribution(type::infected_not_quarantined, simulation_data; output_path="")
    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]
    
    for test_type in keys(simulation_data)
        linestyle = 1
        marker = 1
        labels = Array{String, 1}()
        mytitle = "infected_not_quarantined_dist_$(test_type)_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"
    
        clf()
        figure(figsize=(7,4))
    
        for exp in get!(simulation_data, test_type, Array{Float64, 1}())
            ylim(bottom=0.0, top=1.0)
            plot(exp.second.infected_not_quarantined_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)
    
            xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
            ylabel("Δ = $(exp.second.Δ) hours \n Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)
            title("δ = $(exp.second.δ) minutes", pad=10, fontweight="semibold", fontsize="x-large")
    
            tick_params(labelsize="large")
    
            push!(labels, exp.first)
    
            linestyle = (linestyle + 1) % (length(linestyles)+1)
            marker = (marker + 1) % (length(markers)+1)
    
            if linestyle == 0
                linestyle = 1
            end
            if marker == 0
                marker = 1
            end
        end
        legend(labels, fontsize="large", ncol=2)
        plt.tight_layout()

        savefig("$(output_path)/$(mytitle)")
        println("Infected not quarantined distribution -> saving figure in ... $(output_path)/$(mytitle)\n")
    end    
end


"""
"""
function plot_status_distribution(type::isolation, simulation_data; output_path="")
    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]

    for test_type in keys(simulation_data)
        linestyle = 1
        marker = 1
        labels = Array{String, 1}()
        mytitle = "isolated_status_$(test_type)_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"

        clf()
        figure(figsize=(7,4))

        for exp in get!(simulation_data, test_type, Array{Float64, 1}())
            ylim(bottom=0.0, top=1.0)

            plot(exp.second.isolated_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)
            plot(exp.second.infected_not_isolated_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)
            plot(exp.second.susceptible_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)

            xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
            ylabel("Δ = $(exp.second.Δ) hours \n Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)
            title("δ = $(exp.second.δ) minutes", pad=10, fontweight="semibold", fontsize="x-large")

            tick_params(labelsize="large")

            push!(labels, exp.first)

            linestyle = (linestyle + 1) % (length(linestyles)+1)
            marker = (marker + 1) % (length(markers)+1)

            if linestyle == 0
                linestyle = 1
            end
            if marker == 0
                marker = 1
            end
        end

        labels = ["isolated", "infected not isolated", "susceptibles"]
        legend(labels, fontsize="large", ncol=2)

        plt.tight_layout()
        savefig("$(output_path)/$(mytitle)")
        println("Saving isolated status dist in ... $(output_path)/$(mytitle)\n")
    end
end


"""
"""
function plot_status_distribution(type::quarantine, simulation_data; output_path="")
    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["", "", "", "", "x", "+"]

    for test_type in keys(simulation_data)
        linestyle = 1
        marker = 1
        labels = Array{String, 1}()
        mytitle = "quarantine_status_$(test_type)_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"

        clf()
        figure(figsize=(7,4))

        for exp in get!(simulation_data, test_type, Array{Float64, 1}())
            ylim(bottom=0.0, top=1.0)

            plot(exp.second.quarantined_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)
            plot(exp.second.infected_not_quarantined_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)
            plot(exp.second.susceptible_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)

            xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
            ylabel("Δ = $(exp.second.Δ) hours \n Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)
            title("δ = $(exp.second.δ) minutes", pad=10, fontweight="semibold", fontsize="x-large")

            tick_params(labelsize="large")

            push!(labels, exp.first)

            linestyle = (linestyle + 1) % (length(linestyles)+1)
            marker = (marker + 1) % (length(markers)+1)

            if linestyle == 0
                linestyle = 1
            end
            if marker == 0
                marker = 1
            end
        end

        labels = ["quarantined", "infected not quarantined", "susceptibles"]
        legend(labels, fontsize="large", ncol=2)

        plt.tight_layout()
        savefig("$(output_path)/$(mytitle)")
        println("Saving quarantined status dist in ... $(output_path)/$(mytitle)\n")
    end
end


"""
"""
function plot_SIS_distribution(infected_per_run, susceptible_per_run; output_path="")
    n_intervals = length(infected_per_run[1])

    # interval => [n_infected_iter1...n_infected_itern]
    inf_intervals_to_data = Dict{Int, Array{Float64,1}}()
    sus_intervals_to_data = Dict{Int, Array{Float64,1}}()
    
    for i in 1:n_intervals
        for iter in 1:length(infected_per_run)
            push!(
                get!(inf_intervals_to_data, i, Array{Float64,1}()),
                infected_per_run[iter][i]
            )
    
            push!(
                get!(sus_intervals_to_data, i, Array{Float64,1}()),
                susceptible_per_run[iter][i]
            )
        end
    end
    
    inf_median_to_plot = Array{Float64, 1}()
    inf_25_to_plot = Array{Float64, 1}()
    inf_75_to_plot = Array{Float64, 1}()
    
    sus_median_to_plot = Array{Float64, 1}()
    sus_25_to_plot = Array{Float64, 1}()
    sus_75_to_plot = Array{Float64, 1}()
    
    for i in 1:length(inf_intervals_to_data)
        # infected
        q = quantile(inf_intervals_to_data[i], [0.25, 0.5, 0.75])
    
        push!(inf_median_to_plot, q[2])
        push!(inf_25_to_plot, q[1])
        push!(inf_75_to_plot, q[3])
    
        # susceptible 
        q = quantile(sus_intervals_to_data[i], [0.25, 0.5, 0.75])
    
        push!(sus_median_to_plot, q[2])
        push!(sus_25_to_plot, q[1])
        push!(sus_75_to_plot, q[3])
    end
    
    # PLOT
    clf()
    figure(figsize=(7,4))
    
    #Infected
    plot(inf_median_to_plot, marker="x", markevery=5, markersize=7)
    fill_between(
        range(0, stop=length(inf_median_to_plot)-1), 
        inf_25_to_plot, 
        inf_75_to_plot,
        alpha=0.3
        )
    
    #Suscptible
    plot(sus_median_to_plot, marker="x", markevery=5, markersize=7)
    fill_between(
        range(0, stop=length(sus_median_to_plot)-1), 
        sus_25_to_plot, 
        sus_75_to_plot,
        alpha=0.3
        )
    
    tick_params(labelsize="large")
    ylabel("Quantity", fontweight="roman", fontsize="x-large", labelpad=10)
    xlabel("Time intervals", fontweight="roman", fontsize="x-large", labelpad=10)
    
    legend(["Infected", "Susceptible"])
    
    tight_layout(.5)
    
    my_title = "SIS_$(Dates.format(now(), "Y-mm-ddTHH-MM-SS")).png"
    savefig("$(output_path)/$(my_title)")
    println("SIS -> saving figure in ... $(output_path)/$(my_title)")
end
