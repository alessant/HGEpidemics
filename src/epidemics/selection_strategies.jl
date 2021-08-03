"""
    uniform(h::Hypergraph, α::Float64; kwargs...)

 Select a random sample of `α` nodes or `α` hyperdeges of the
 hypergraph `h` to immunize.
"""
function uniform(h::Hypergraph, α::Union{Int, Float64}; kwargs...)
    rng = MersenneTwister(1234)
    to_return = shuffle!(rng, collect(1:nhv(h)))
    to_return[1:ceil(Int, length(to_return)*α)]
end

"""
    Return top-α elements pre-computed.
"""
function custom(h::Hypergraph, α::Union{Int, Float64, Nothing}; kwargs...)
    if !isnothing(α)
        to_return = deserialize(kwargs[:path])
        return to_return[1:ceil(Int, length(to_return)*α)]
    end
    return deserialize(kwargs[:path])
end