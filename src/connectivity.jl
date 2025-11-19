#=
Module for creating adjacency matrices 
=#

using Geodesy, SparseArrays
using Base.Threads


function balloon_distance(pos1, pos2)  # pos = [lat, lon, alt]
    lla1 = LLA(pos1[1], pos1[2], pos1[3] * 1000)  # altitude in meters
    lla2 = LLA(pos2[1], pos2[2], pos2[3] * 1000)
    
    ecef1 = ECEF(lla1, wgs84)
    ecef2 = ECEF(lla2, wgs84) 
    
    return euclidean_distance(ecef1, ecef2) / 1000  # Convert to km
end


# Computes weighting for connectivity graph
function weight(dist_km, min_cutoff=100, max_cutoff=300)::Float64
    if dist_km < min_cutoff
        return 1.0
    else
        return max(0.0,(max_cutoff - dist_km)/(max_cutoff - min_cutoff))
    end
end



function build_adjacency!(network::Network)
    n = network.n_balloons
    
    # Pre-allocate arrays for each thread to avoid race conditions
    n_threads = Threads.nthreads()
    I_arrays = [Int[] for _ in 1:n_threads]
    J_arrays = [Int[] for _ in 1:n_threads] 
    V_arrays = [Float64[] for _ in 1:n_threads]
    
    # Parallel loop over i
    Threads.@threads for i in 1:n
        thread_id = Threads.threadid()
        I_local = I_arrays[thread_id]
        J_local = J_arrays[thread_id]
        V_local = V_arrays[thread_id]
        
        for j in (i+1):n
            pos_i = network.positions[i, :]
            pos_j = network.positions[j, :]
            
            dist = balloon_distance(pos_i, pos_j)
            w = weight(dist)
            
            if w > 0.0
                push!(I_local, i); push!(J_local, j); push!(V_local, w)
                push!(I_local, j); push!(J_local, i); push!(V_local, w)
            end
        end
    end
    
    # Combine results from all threads
    I = vcat(I_arrays...)
    J = vcat(J_arrays...)
    V = vcat(V_arrays...)
    
    network.adjacency = sparse(I, J, V, n, n)
     build_laplacian_sparse!(network)
end