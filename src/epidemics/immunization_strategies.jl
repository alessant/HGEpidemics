# """
#     uniform(h::Hypergraph, α::Float64; kwargs...)

#  Select a random sample of `α` nodes or `α` hyperdeges of the
#  hypergraph `h` to immunize.
# """
# function uniform(h::Hypergraph, α::Union{Int, Float64}; kwargs...)
#     rng = MersenneTwister(1234)
#     to_return = shuffle!(rng, collect(1:nhv(h)))
#     to_return[1:ceil(Int, length(to_return)*α)]
# end



"""
    degrees(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to their degree in `h`.

 The algorithm selects `α` nodes with the higher degree.
"""
function degrees(h::Hypergraph, α::Float64; kwargs...)
    #d = [length(gethyperedges(h, v)) for v in 1:nhv(h)]
    g, gw, hg, hgw = generate_hg_g_weighted!(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:user2vertex],
        kwargs[:loc2he],
        kwargs[:δ]
    )

    d = [length(gethyperedges(hg, v)) for v in 1:nhv(hg)]

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



"""
    centrality(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to their betweeness centrality in `h`.

 The algorithm selects `α` nodes with the higher bc values.
"""
function centrality(h::Hypergraph, α::Float64; kwargs...)
    m = SimpleHypergraphs.adjacency_matrix(h; s=2)

    g = LightGraphs.SimpleGraph(m)
    bc = LightGraphs.betweenness_centrality(g)

    to_return = sort!(collect(1:nhv(h)), by = x -> bc[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



"""
    random_walk(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to their random walk centrality in `h`.

 The algorithm selects `α` nodes with the higher rw values.
"""
function random_walk(h::Hypergraph, α::Float64; kwargs...)
    rwc = Dict{Int, Int}(v => 0 for v in 1:nhv(h))

    for r = 1:100
        for v = 1:nhv(h)
            length(gethyperedges(h, v)) == 0 && continue
            rwc[SimpleHypergraphs.random_walk(h, v)] += 1
        end
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> rwc[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



"""
    acquaintance(h::Hypergraph, α::Float64; kwargs...)

 Select `α` nodes or `α` hyperdeges of the hypergraph `h` to immunize
 according to the acquaintance strategy. It also makes use of local
 knowledge selecting the neighbor with higher degree.
"""
function acquaintance(h::Hypergraph, α::Float64; kwargs...)
    rng = MersenneTwister(1234)
    nodes = shuffle!(rng, collect(1:nhv(h)))

    n_to_immunize = ceil(Int, length(nodes)*α)
    to_immunize = Set{Int}()

    for n in nodes
        neighbors = Set{Int}()
        for he in keys(h.v2he[n])
            union!(neighbors, keys(h.he2v[he]))
        end
        delete!(neighbors, n) #remove v from its neighborhood

        #consider only the nodes that have at least one neighbor
        length(neighbors) == 0 && continue

        # we want to select the neighbor with higher degree
        n_degrees = [(length(gethyperedges(h, n)), n) for n in neighbors]
        chosen = maximum(n_degrees)[2]

        push!(to_immunize, chosen)

        length(to_immunize) == n_to_immunize && break
    end

    return collect(to_immunize)
end



# """
#     lockdown(h::Hypergraph, α::Union{Int, Float64}; kwargs...)

#  Close all location indicated in `kwargs[:path]`.
# """
# function lockdown(h::Hypergraph, α::Union{Int, Float64, Nothing}; kwargs...)
#     if !isnothing(α)
#         to_return = deserialize(kwargs[:path])
#         return to_return[1:ceil(Int, length(to_return)*α)]
#     end
#     return deserialize(kwargs[:path])
# end



function hg_pagerank(h::Hypergraph, α::Union{Int, Float64, Nothing}; kwargs...)
    g, gw, hg, hgw = generate_hg_g_weighted!(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:user2vertex],
            kwargs[:loc2he],
            kwargs[:δ]
        )

    for v in 1:nhv(hg)
        for he in 1:nhe(hg)
            if isnothing(hg[v, he])
                hg[v, he] = 0
            end
        end
    end

    R, P = hrwr(hg)
    rank = R[1, :]

    to_return = sort!(collect(1:nhv(h)), by = x -> rank[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end


function whg_pagerank(h::Hypergraph, α::Union{Int, Float64, Nothing}; kwargs...)
    g, gw, hg, hgw = generate_hg_g_weighted!(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:user2vertex],
            kwargs[:loc2he],
            kwargs[:δ]
        )

    R, P = hrwr(hgw)
    rank = R[1, :]

    to_return = sort!(collect(1:nhv(h)), by = x -> rank[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end


function g_pagerank(h::Hypergraph, α::Union{Int, Float64, Nothing}; kwargs...)
    g, gw, hg, hgw = generate_hg_g_weighted!(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:user2vertex],
            kwargs[:loc2he],
            kwargs[:δ]
        )

    rank = pagerank(g)

    to_return = sort!(collect(1:nhv(h)), by = x -> rank[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end


function wg_pagerank(h::Hypergraph, α::Union{Int, Float64, Nothing}; kwargs...)
    g, gw, hg, hgw = generate_hg_g_weighted!(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:user2vertex],
            kwargs[:loc2he],
            kwargs[:δ]
        )

    rank = weighted_pagerank(gw)

    to_return = sort!(collect(1:nhv(h)), by = x -> rank[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end


# entropy
function hg_entropy(h::Hypergraph, α::Float64; kwargs...)
    # g, gw, hg, hgw = generate_hg_g_weighted!(
    #     kwargs[:df],
    #     kwargs[:mindate],
    #     kwargs[:maxdate],
    #     kwargs[:user2vertex],
    #     kwargs[:loc2he],
    #     kwargs[:δ];
    #     all = kwargs[:history]
    # )

    # dict
    entropies = evaluate_entropy(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:user2vertex],
        kwargs[:loc2he]
        )

    d = Array{Float64, 1}()

    for v in 1:nhv(h)
        push!(d, entropies[v])
    end

    # for v in 1:nhv(hgw)
    #     direct = 0.0
    #     for he in gethyperedges(hgw, v)
    #         direct += hgw[v, he.first]
    #     end
    #     direct /= length(gethyperedges(hgw, v))

    #     #val = direct * entropies[v]
    #     val = entropies[v]
    #     push!(d, val)
    # end

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



function hg_entropy_mean(h::Hypergraph, α::Float64; kwargs...)
    # dict
    entropies = evaluate_entropy_indirect_mean(
        kwargs[:df],
        kwargs[:intervals],
        kwargs[:user2vertex],
        kwargs[:loc2he],
        kwargs[:start_slot],
        kwargs[:end_slot]
        )

    d = Array{Float64, 1}()

    for v in 1:nhv(h)
        push!(d, entropies[v])
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



function g_degrees(h::Hypergraph, α::Float64; kwargs...)
    g, gw, hg, hgw = generate_hg_g_weighted!(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:user2vertex],
        kwargs[:loc2he],
        kwargs[:δ]
    )

    d = degree(g)

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



function g_entropy(h::Hypergraph, α::Float64; kwargs...)

    # dict
    entropies = evaluate_entropy_direct(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:δ],
        kwargs[:user2vertex],
        kwargs[:loc2he]
        )

    d = Array{Float64, 1}()

    for v=1:nhv(h)
        if haskey(entropies, v)
            push!(d, entropies[v])
        else
            push!(d, typemin(Int))
        end
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end




function g_entropy_mean(h::Hypergraph, α::Float64; kwargs...)
    # dict
    entropies = evaluate_entropy_direct_mean(
        kwargs[:df],
        kwargs[:intervals],
        kwargs[:δ],
        kwargs[:user2vertex],
        kwargs[:loc2he],
        kwargs[:start_slot],
        kwargs[:end_slot]
        )

    d = Array{Float64, 1}()

    for v=1:nhv(h)
        if haskey(entropies, v)
            push!(d, entropies[v])
        else
            push!(d, typemin(Int))
        end
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end



function both_entropy(h::Hypergraph, α::Float64; kwargs...)
    # dict
    entropies_direct = evaluate_entropy_direct(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:δ],
        kwargs[:user2vertex],
        kwargs[:loc2he]
    )

    entropies_indirect = evaluate_entropy(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:user2vertex],
        kwargs[:loc2he]
    )

    d = Array{Float64, 1}()

    for v=1:nhv(h)
        if haskey(entropies_direct, v)
            push!(d, (entropies_direct[v]) * (entropies_indirect[v]))
        else
            push!(d, eps() * (entropies_indirect[v]))
        end
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end




function both_entropy_norm(h::Hypergraph, α::Float64; kwargs...)
        #println(slots)
    # 
    # dict
    entropies_direct = evaluate_entropy_direct(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:δ],
        kwargs[:user2vertex],
        kwargs[:loc2he]
    )

    entropies_indirect = evaluate_entropy(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:user2vertex],
        kwargs[:loc2he]
    )

    d = Array{Float64, 1}()

    mean_direct = mean(values(entropies_direct))
    mean_indirect = mean(values(entropies_indirect))

    for v=1:nhv(h)
        if haskey(entropies_direct, v)
            push!(d, (entropies_direct[v] / mean_direct) + (entropies_indirect[v] / mean_indirect))
        else
            push!(d, 0 + (entropies_indirect[v] / mean_indirect))
        end
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end




function both_entropy_mean(h::Hypergraph, α::Float64; kwargs...)
    # dict
    entropies_direct = evaluate_entropy_direct_mean(
        kwargs[:df],
        kwargs[:intervals],
        kwargs[:δ],
        kwargs[:user2vertex],
        kwargs[:loc2he],
        kwargs[:start_slot],
        kwargs[:end_slot]
    )

    entropies_indirect = evaluate_entropy_indirect_mean(
        kwargs[:df],
        kwargs[:intervals],
        kwargs[:user2vertex],
        kwargs[:loc2he],
        kwargs[:start_slot],
        kwargs[:end_slot]
    )

    d = Array{Float64, 1}()

    for v=1:nhv(h)
        if haskey(entropies_direct, v)
            push!(d, entropies_direct[v] * entropies_indirect[v])
        else
            push!(d, eps() * entropies_indirect[v])
        end
    end

    to_return = sort!(collect(1:nhv(h)), by = x -> d[x], rev = true)
    to_return[1:ceil(Int, length(to_return)*α)]
end






function acquaintance_new(h::Hypergraph, α::Float64; kwargs...)
    rng = MersenneTwister(1234)
    nodes = shuffle!(rng, collect(1:nhv(h)))

    n_to_immunize = ceil(Int, length(nodes)*α)
    to_immunize = Set{Int}()

    g = nothing
    entropies = nothing

    g, gw, hg, hgw = generate_hg_g_weighted!(
        kwargs[:df],
        kwargs[:mindate],
        kwargs[:maxdate],
        kwargs[:user2vertex],
        kwargs[:loc2he],
        kwargs[:δ]
    )

    if kwargs[:type] == :entropy
        entropies = evaluate_entropy(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:user2vertex],
            kwargs[:loc2he]
            )
    elseif kwargs[:type] == :entropy_g
        entropies = evaluate_entropy_direct(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:δ],
            kwargs[:user2vertex],
            kwargs[:loc2he]
            )
    elseif kwargs[:type] == :entropy_both
        entropies_direct = evaluate_entropy_direct(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:δ],
            kwargs[:user2vertex],
            kwargs[:loc2he]
        )
    
        entropies_indirect = evaluate_entropy(
            kwargs[:df],
            kwargs[:mindate],
            kwargs[:maxdate],
            kwargs[:user2vertex],
            kwargs[:loc2he]
        )
    
        entropies = Array{Float64, 1}()
    
        for v=1:nhv(h)
            if haskey(entropies_direct, v)
                push!(entropies, entropies_direct[v] * entropies_indirect[v])
            else
                push!(entropies, eps() * entropies_indirect[v])
            end
        end
    else
        #nothing
    end

    for n in nodes
        chosen = nothing

        # neighbors = Set{Int}()
        # for he in keys(hg.v2he[n])
        #     union!(neighbors, keys(hg.he2v[he]))
        # end
        # delete!(neighbors, n) #remove v from its neighborhood
        neighbors = LightGraphs.neighbors(g, n)

        #consider only the nodes that have at least one neighbor
        length(neighbors) == 0 && continue

        # we want to select the neighbor 
        # according to a given strategy
        if kwargs[:type] == :uniform
            neighbors = collect(neighbors)
            chosen = neighbors[rand(1:length(neighbors))]

        elseif kwargs[:type] == :degree
            #with higher degree
            n_degrees = [(length(gethyperedges(h, n)), n) for n in neighbors]
            chosen = maximum(n_degrees)[2]
        
        elseif kwargs[:type] == :degree_g
            n_degrees = [(length(LightGraphs.neighbors(g, n)), n) for n in neighbors]
            chosen = maximum(n_degrees)[2]

        elseif kwargs[:type] == :entropy
            n_entropies = [(entropies[n], n) for n in neighbors]
            chosen = maximum(n_entropies)[2]

        elseif kwargs[:type] == :entropy_g
            n_entropies = [(entropies[n], n) for n in neighbors]
            chosen = maximum(n_entropies)[2]
        else #entropy both
            n_entropies = [(entropies[n], n) for n in neighbors]
            chosen = maximum(n_entropies)[2]
        end

        push!(to_immunize, chosen)
        length(to_immunize) == n_to_immunize && break
    end

    return collect(to_immunize)
end



function most_loved(h::Hypergraph, α::Union{Int, Float64, Nothing}; kwargs...)
    return deserialize(kwargs[:path])
end