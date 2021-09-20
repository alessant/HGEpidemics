function brwr(h::Hypergraph, heweights = nothing)
    H = fill(0, nhv(h), nhe(h))
    Dᵥ, W = get_weighted_v_degree_matrix!(h, H)

    #users x users
    A₁ = fill(0, nhe(h), nhe(h))
    #users x items
    A₂ = transpose(W)
    #items x users
    A₃ = W
    #items x items
    A₄ = fill(0, nhv(h), nhv(h))

    A = vcat(hcat(A₁, A₂), hcat(A₃, A₄))
    D = initializeD(A)

    P = transpose(inv(D) * A)

    S = transpose(W)#computeSimilarity(transpose(W))


    c = .8
    ϵ = .0001
    nelem = nhv(h) + nhe(h)

    user_to_items_rank = Dict{Int, Any}()

    for user=1:nhe(h)
        println("User $(user)")

        q = fill(1/nelem, nelem)
        #q = fill(0.0, nelem)
        #q[user] = 1


        # for v=1:nhe(h)
        #     w = 0.0
        #     for e in gethyperedges(h,v)
        #         w += S[user, e.first] * h[v, e.first]
        #     end
        #     q[nhv(h) + v] = w
        # end
        #
        # q = q ./ sum(q)

        r = fill(1/nelem, nelem)
        r_prev = fill(1/nelem, nelem)

        while true
            for i=1:nelem
                r[i] = (c * dot(P[:,i], r_prev)) + (1-c) * q[i]
            end

            if norm(r_prev - r) < ϵ
                break
            end

            r_prev = copy(r)
        end

        push!(
            user_to_items_rank,
            user => r[nhe(h)+1 : nelem]
        )

    end

    #H, Dᵥ, Dₑ, Wₑ, W, P, user_to_items_rank
    A, W, D, user_to_items_rank
end


function bevaluate(h, ranking_matrix)
    u_metrics = Dict{Int, Any}()

    #for each user
    for he=1:nhe(h)
        seed_nodes = getvertices(h, he)
        u_ranks = get!(ranking_matrix, he, nothing) #vector
        print(seed_nodes)
        obs = fill(0.0, length(seed_nodes))
        i = 1
        for v in seed_nodes
            obs[i] = h[v.first, he]
            i=i+1
        end

        j=1

        prs = fill(0.0, length(seed_nodes))
        #for each item of u
        i = 1
        for u in seed_nodes
            prs[i] = u_ranks[u.first]
            i=i+1
        end

        #add corr to corrs
        rho = corspearman(obs, prs)

        j+=1

        push!(
            u_metrics,
            he => [rho]
        )

    end

    u_metrics
end


function initializeD(A)
    D = fill(0.0, size(A))

    for row=1:size(A)[1]
        D[row, row] = sum(A[row, :])
    end

    D
end


function get_weighted_v_degree_matrix!(h, H, heweights = nothing)
    D = fill(0.0, nhv(h), nhv(h))
    W = fill(0.0, nhv(h), nhe(h))

    if isnothing(heweights)
        heweights = fill(1, nhe(h))
    end

    for v=1:nhv(h)
        wᵥ = 0.0
        for he in gethyperedges(h, v)
            wᵥ += heweights[he.first]
            H[v, he.first] = 1
            W[v, he.first] = h[v, he.first]
        end
        D[v, v] = wᵥ
    end

    D, W
end


function get_weighted_he_degree_matrix(h, heweights = nothing)
    if isnothing(heweights)
        heweights = fill(1, nhe(h))
    end

    D = fill(0.0, nhe(h), nhe(h))
    W = fill(0.0, nhe(h), nhe(h))

    for he=1:nhe(h)
        wₑ = 0.0
        for v in getvertices(h, he)
            wₑ += h[v.first, he]
        end

        D[he, he] = wₑ #weights(he) * length(getvertices(h, he))
        W[he, he] = heweights[he]
    end

    D, W
end
