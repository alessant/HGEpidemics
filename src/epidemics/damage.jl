"""
    check damage
"""

function store_damage(
    sim_data::Array, 
    user2vertex::Dict, 
    output_path::AbstractString
)   

    vertex2user = Dict{Int, String}(values(user2vertex) .=> keys(user2vertex))

    open(output_path, "w+") do file
        write(file, "userid,vertexid,value\n")

        for i=1:length(sim_data)
            write( 
                file,
                "$(vertex2user[i]),$i,$(sim_data[i])\n"
            )
        end
    end
end


function check_damage(
    sim_data::Array,
    ground_truth_path::AbstractString
)

    gt_vals = CSV.read(ground_truth_path, DataFrame).value
    tmp = [(gt_vals[i]-sim_data[i])/gt_vals[i] for i in 1:length(sim_data) if gt_vals[i] > 0]
    sum(tmp)/length(tmp)
end
