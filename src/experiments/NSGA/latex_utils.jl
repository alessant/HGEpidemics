function save_spider_data(output_path, min_infected_conf, min_damage_conf, min_distance_conf)

    names = "Features\tmin_infected\tmin_damage\tmin_distance\n"

    function to_boolean(val)
        println(val)
        val >= 0.5 ? 1.0 : 0.0
    end

    open(output_path, "w+") do f
        write(f, names)

        write(f, "PPM\t$(min_infected_conf[1]*10)\t$(min_damage_conf[1]*10)\t$(min_distance_conf[1]*10)\n")

        write(
            f, 
            "EM\t$(to_boolean(min_infected_conf[2])*10)\t$(to_boolean(min_damage_conf[2])*10)\t$(to_boolean(min_distance_conf[2])*10)\n"
        )
        
        write(f, "Quarantine\t$(min_infected_conf[4]*10)\t$(min_damage_conf[4]*10)\t$(min_distance_conf[4]*10)\n")
        write(f, "Tracing\t$(min_infected_conf[6]*10)\t$(min_damage_conf[6]*10)\t$(min_distance_conf[6]*10)\n")
        write(f, "Isolation\t$(min_infected_conf[3]*10)\t$(min_damage_conf[3]*10)\t$(min_distance_conf[3]*10)\n")

        write(
            f, 
            "{Avoiding Crowding}\t$(to_boolean(min_infected_conf[7])*10)\t$(to_boolean(min_damage_conf[7])*10)\t$(to_boolean(min_distance_conf[7])*10)\n"
        )

        write(f, "Lockdown\t$(min_infected_conf[5]*10)\t$(min_damage_conf[5]*10)\t$(min_distance_conf[5]*10)\n")
    end
end


function save_scatter_data(output_path, population, colors)
    header = "x,y,class\n"

    open(output_path, "w+") do f
        write(f, header)

        for index in 1:length(population)
            to_write = "$(population[index].y[1]),$(population[index].y[2]),$(colors[index])\n"
            write(f, to_write)
        end
    end
end


function save_3Dscatter_data(output_path, population, colors)
    header = "x,y,z,class\n"

    open(output_path, "w+") do f
        write(f, header)

        for index in 1:length(population)
            to_write = "$(population[index].y[1]),$(population[index].y[2]),$(population[index].y[3]),$(colors[index])\n"
            write(f, to_write)
        end
    end
end
