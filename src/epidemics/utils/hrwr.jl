"""
    hrwr(h::Hypergraph; heweights = nothing)

Compute a rank matrix `R`, where each cell `(i, j)` is the
rank value of the node `j` with respect to node `i`,
obtained with a random walk with restart on the hypergraph `h`.

For more information, see the following paper:
Abdelghani Bellaachia and Mohammed Al-Dhelaan,
_Random Walks in Hypergraphs_
Proceedings of the 2013 International Conference on
Applied Mathematics and Computational Methods.

# Optional arguments:
- `heweights`: hyperedge weights
- `α`: dumping factor
- `ϵ`: convergence threshold
- `iter`: number of iterations
"""
function hrwr(h::Hypergraph; heweights=nothing, α =.85, ϵ = 1.0e-6, iter=100)
    H, W, Dᵥ = compute_v_matrices!(h; heweights=heweights)
    Wₑ, Dₑ = compute_he_matrices(h; heweights=heweights)

    P = inv(Dᵥ) * H
    P = P * Wₑ
    P = P * inv(Dₑ)
    P = P * transpose(W)

    Pₜ = transpose(P)

    R = fill(0.0, nhv(h), nhv(h))

    for i=1:nhv(h)
        q = fill(1/nhv(h), nhv(h))
        #q = fill(0, nhv(h))
        #q[i] = 1
        r = fill(1/nhv(h), nhv(h))
        r_prev = fill(1/nhv(h), nhv(h))

        liter = iter
        while liter >= 0
            # for i=1:nhv(h)
            #     r[i] = (α * dot(P[:, i], r_prev)) + (1-α) * q[i]
            # end
            r = α * (Pₜ * r_prev) + (1-α) * q

            if norm(r_prev - r) < ϵ
                break
            end
            r_prev = deepcopy(r)
            liter -= 1
        end
        if liter == -1
            error("RWR did not converge after $iter iterations.")
        end
        R[i, :] = copy(r)
    end

    R, P
end


"""
    compute_v_matrices!(h; heweights = nothing)

Auxiliary function for random walk with restart on hypergraphs.
Computes:
    - `H`: the incidence matrix of the hypergraph `h`
    - `W`: the weighted incidence matrix of the hypergraph `h`
    - `D`: the diagonal matrix of the weighted degree of vertices
"""
function compute_v_matrices!(h; heweights=nothing)
    H = fill(0, nhv(h), nhe(h))
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

        if wᵥ == 0
            D[v, v] = eps()
        else
            D[v, v] = wᵥ
        end
    end
    H, W, D
end


"""
    compute_he_matrices(h, heweights = nothing)

Auxiliary function for random walk with restart on hypergraphs.
    Computes:
    - `Wₑ`: diagonal matrix of hyperedge weights
    - `Dₑ`: diagonal matrix for weighted degree of a hyperedge
"""
function compute_he_matrices(h; heweights=nothing)
    D = fill(0.0, nhe(h), nhe(h))
    W = fill(0.0, nhe(h), nhe(h))

    if isnothing(heweights)
        heweights = fill(1, nhe(h))
    end

    for he=1:nhe(h)
        wₑ = 0.0
        for v in getvertices(h, he)
            wₑ += h[v.first, he]
        end

        if wₑ == 0
            D[he, he] = eps()
        else
            D[he, he] = wₑ
        end
        W[he, he] = heweights[he]
    end

    W, D
end
