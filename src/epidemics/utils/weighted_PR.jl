function weighted_pagerank(
    g::AbstractGraph{U}, 
    α=0.85, 
    n::Integer=100,
    ϵ=1.0e-6
    ) where U <: Integer

    α_div_outdegree = Vector{Float64}(undef, nv(g))
    den_weighted_degree = Vector{Float64}(undef, nv(g))
    dangling_nodes = Vector{U}()

    for v in vertices(g)
        if outdegree(g, v) == 0
            push!(dangling_nodes, v)
        end
        α_div_outdegree[v] = (α/outdegree(g, v))

        w = 0.0

        for u in neighbors(g, v)
            w += g.weights[u, v]
        end

        den_weighted_degree[v] = w
    end

    N = Int(nv(g))
    # solution vector and temporary vector
    x = fill(1.0 / N, N)
    xlast = copy(x)
    for _ in 1:n
        dangling_sum = 0.0
        for v in dangling_nodes
            dangling_sum += x[v]
        end
        # flow from teleprotation
        for v in vertices(g)
            xlast[v] = (1 - α + α * dangling_sum) * (1.0 / N)
        end

        # flow from edges
        for v in vertices(g)
            for u in inneighbors(g, v)
                #println("$v -- $u -- $(g.weights[u, v]) / $(den_weighted_degree[u])")
                xlast[v] += (x[u] * α * (g.weights[u, v] / den_weighted_degree[u])) #α_div_outdegree[u])
            end
        end
        # l1 change in solution convergence criterion
        err = 0.0
        for v in vertices(g)
            err += abs(xlast[v] - x[v])
            x[v] = xlast[v]
        end
        if (err < N * ϵ)
            return x
        end
    end
    error("Pagerank did not converge after $n iterations.") 
end