####################################################
# ALGORITHM TO GENERATE WEIGHTS FOR TANK INVENTORIES 
#################################################### 

## STRUCTURE
# P - Identify Flow from Producer Nodes
# C - Identify Flows to Consumer Nodes 
# alpha_f - Alpha for Outflows to Consumer Nodes (To Calculate Delta C)
# EI - Edge Weights (Weights assigned to each Flow Rates (To Calculate Distance from Storage to Consumer) 
# NI - Number of Inventories 
# O - Number of Flows 
# M - Connectivity Matrix 

function assign_weights()

    ##############################
    # DATA FOR WEIGHT ASSIGNMENT  
    ############################## 
    # nc - number of consumer nodes 
    # np - number of producer nodes 
    # M_algo - Matrix includig producer and consumer nods 
    # CINDEX - Indexes for consumer nodes in M algo 
    np = size(P)[1]
    nc = size(C)[1]
    M_algo = vcat(M, P)
    M_algo = vcat(M_algo, C) 
    CINDEX  = NI + np+1:NI + np + nc

    ##############################
    # CALCULATE DELTA C
    ############################## 
    # nc_idx - indexes for consumer nodes 
    # deltac - for each edge connected to a consumer node, calculate deltac 
    nc_idx = zeros(nc)
    for k in 1:nc
        idx = findfirst(x -> x == 1, C[k, :])
        nc_idx[k] = idx
    end
    nc_idx = Int.(nc_idx)

    deltac = zeros(nc) 
    for x in 1:nc
        deltac[x] = 1/alpha_f[nc_idx[x]]
    end

    ##############################
    # CREATE WEIGHTED DIGRAPH
    ############################## 
    # n,m - Nodes (n) and Edges (m) in the Network
    # G - Directed graph 

    n, m = size(M_algo)
    G = SimpleWeightedDiGraph(n)

    # Build network from incidence matrix
    for j in 1:m # Go through each column in the M_algo matrix 
        start_node = findfirst(x -> x == -1, M_algo[:, j])
        end_node   = findfirst(x -> x == 1, M_algo[:, j])
        if !isnothing(start_node) && !isnothing(end_node)
            add_edge!(G, start_node, end_node, EI[j])
        end
    end

    ##############################
    # VISUALIZE NETWORK 
    ############################## 
    grph = gplot(G)
    display(grph)

    ##############################
    # CALCULATE DCS
    ############################## 
    # dcs - shortest distance from each inventory to the consumer node 
    dcs = zeros(nc, NI)
    for x in 1:nc
        for y in 1:NI
            distance = dijkstra_shortest_paths(G, y)
            dcs[x, y] = distance.dists[CINDEX[x]]
        end
    end  
    alpha_inventory = [1/minimum(col) for col in eachcol(dcs.+deltac)] 

    return alpha_inventory

end



