using JSON
using Dates
using Serialization

input = "data/got/episodes.json"
output = "data/got/checkins_JCS.csv"

save_checkins(input, output)


function save_checkins(input, output)
    gotjson = JSON.parsefile(input; dicttype=Dict, inttype=Int64, use_mmap=true)

    clock = Dates.DateTime(2020,01,01,00,00,00) #1577833200000

    # chars = deserialize("src/experiments/immunization/got/data/top_characters.data")
    # println(chars)

    open(output, "w+") do file
        header = "character,location,timestamp"
        println(file, header)
        
        day = 0

        for i in range(1, stop=11)

            for episode in gotjson["episodes"] 
                last_scene = nothing
                
                for scene in episode["scenes"]
                    
                    #locations
                    location = nothing

                    if haskey(scene, "subLocation")
                        location = scene["subLocation"]
                    else
                        location = scene["location"]
                    end

                    sceneStart = Dates.Time(scene["sceneStart"], "H:M:S")
                    sceneEnd = Dates.Time(scene["sceneEnd"], "H:M:S")
                    if isnothing(last_scene)
                        #timestamp
                        timestamp = clock + Hour(sceneEnd) + Minute(sceneEnd) + Second(sceneEnd)
                        clock = timestamp

                    else
                        prevStart =  Dates.Time(last_scene["sceneStart"], "H:M:S")

                        #elapsed =  sceneStart
                        #          - Hour(prevStart) - Minute(prevStart) - Second(prevStart)
                        elapsed = sceneEnd - Hour(sceneStart) -  Minute(sceneStart) - Second(sceneStart)
                                
                        timestamp = clock + Hour(elapsed) + Minute(elapsed) + Second(elapsed)
                        clock = timestamp
                    end 
                
                    
                    #characters
                    for c in scene["characters"]
                        name = Symbol(replace(c["name"], " " => "_"))

                        # !(name in chars) && continue

                        to_write = string(name, ",", location, ",", timestamp)

                        println(file, to_write)
                    end
                    last_scene = scene
                end
            end     
        end #end repetition
    end
end



# sceneStart = Dates.Time("0:00:40", "H:M:S")
# elapsed =  Dates.Time("0:01:45", "H:M:S")

# elapsed - Hour(sceneStart) - Minute(sceneStart) - Second(sceneStart)


# d2 = Dates.DateTime(2020,01,02,00,10,00)

# r = d2 - basestamp
# typeof(r)

# convert(DateTime, r)


# basestamp = Dates.DateTime(2020,01,01,00,00,00)
# t1 = Dates.Time("00:10:00", "H:M:S")
# t2 = Dates.Time("00:00:10", "H:M:S")

# step1 = Dates.value(basestamp) +  Dates.value(t1)

# basestamp = convert(DateTime, Millisecond(step1))

# step2 = Dates.value(basestamp) + Dates.value(t1)

# convert(DateTime, Millisecond(step2))


# Dates.format(convert(DateTime, period), "MM:SS")

# t1 = Dates.Time("0:00:40", "H:M:S")
# t2 = Dates.Time("0:01:45", "H:M:S")

# Dates.value(t1)

# res = DateTime(Date(basestamp), t1)
# res2 = DateTime(Date(res), t2)



# basestamp = Dates.DateTime(2020,01,01,00,00,00)
# t1 = Dates.Time("00:10:00", "H:M:S")
# h, m, s = Hour(t1), Minute(t1), Second(t1)

# basestamp = basestamp + h + m + s
# basestamp = basestamp + h + m + s

# for i=1:100
#     global basestamp
#     basestamp = basestamp + h + m + s
# end
# basestamp

# Dates.value(basestamp) + 
# Millisecond(t1)

# convert(DateTime, Millisecond(Dates.value(t1)))

# convert(Millisecond, Dates.value(t1))

# Dates.format(convert(DateTime, Dates.value(t1)), "MM:SS")