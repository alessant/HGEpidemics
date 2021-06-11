"""
    Utilities functions to store the output 
    in a format useful to plot data in Latex
"""
function store_infected_distribution_data(simulation_data; output_path="")
       
    for k in keys(simulation_data)
        to_write = Dict{String, Array{Float64}}() 
        exp_data = simulation_data[k]

        for elem in exp_data
            label = elem.first
            data = elem.second.infected_distribution
            
            push!(
                to_write,
                label => data
            )
        end

        my_title = "$(output_path)/latex-infected-dist.csv"
        println("Latex infected distribution -> saving data in ... $(my_title)")

        open(my_title, "w+") do f
            #header
            write(f, "x,")
            header = collect(keys(to_write))

            for (index,elem) in enumerate(header)
                index == length(header) ? write(f, "$elem\n") : write(f, "$elem,")
            end

            #core
            for i in 1:length(collect(values(to_write))[1])
                write(f, "$i,")
                for (index,elem) in enumerate(header)
                    index == length(header) ? write(f, "$(to_write[elem][i])\n") : write(f, "$(to_write[elem][i]),")
                end
            end
        end 
    end
end


"""
"""
function store_SIS_distribution_data(infected_per_run, susceptible_per_run; output_path="")
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

    my_title = "$(output_path)/latex-SIS-dist.csv"
    println("Latex infected distribution -> saving figure in ... $(my_title)")

    open(my_title, "w+") do f
        write(f, "x,infected,infected-25,infected-75,susceptible,susceptible-25,susceptible-75\n")
    
        for i in 1:length(inf_median_to_plot)
            infected = inf_median_to_plot[i]
            infected_25 = inf_25_to_plot[i]
            infected_75 = inf_75_to_plot[i]
    
            susceptible = sus_median_to_plot[i]
            susceptible_25 = sus_25_to_plot[i]
            susceptible_75 = sus_75_to_plot[i]
    
            write(
                f,
                "$i,$infected,$infected_25,$infected_75,$susceptible,$susceptible_25,$susceptible_75\n"
            )
        end
    end
end
