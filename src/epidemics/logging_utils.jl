"""
 Simulation logging utils.
"""
function init_log(path::AbstractString, header::AbstractString)
    open(path, "w+") do fhandle
        write(fhandle, header)
    end
end


function log(path::AbstractString, to_write::AbstractString)
    open(path, "a") do fhandle
        write(fhandle, to_write)
    end
end