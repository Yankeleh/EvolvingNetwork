#=
Network struct, constructor, simulate!, helpers
=#


using LinearAlgebra, SparseArrays, Random
using StatsBase, Printf


#Structure for keeping track of spectral properties for network

mutable struct SpectralHistory
    algebraic_connectivity::Vector{Float64}
    num_components::Vector{Int}
    largest_component_size::Vector{Int}
    cheeger_lower_bound::Vector{Float64}
    cheeger_upper_bound::Vector{Float64}
    timestamps::Vector{Float64}
end


#Network structure

mutable struct Network
    net_size::Int
    positions::Matrix{Float64} #net_size x 3 (lat_deg, lon_deg, alt_km)
    
    adjacency::SparseMatrixCSC{Float64, Int}
    laplacian::SparseMatrixCSC{Float64, Int}
    
    spectral_history::SpectralHistory
end


#Constructor for the network


function Network(net_size::Int)
    # Initialize positions using proper sphere sampling

    positions = Matrix{Float64}(undef, net_size, 3)

    for i in 1:net_size
        # Uniform sampling on sphere surface
        θ = 2π * rand()  # Longitude: 0 to 2π
        φ = acos(1 - 2*rand())  # Latitude: compensates for sphere curvature
        
        # Convert to degrees and set altitude
        lat_deg = (π/2 - φ) * 180/π  
        lon_deg = θ * 180/π - 180    
        alt_km = 10.0
        
        positions[i, :] = [lat_deg, lon_deg, alt_km]
    end
    
    # Initialize empty matrices (will be built by build_adjacency!)
    adjacency = spzeros(Float64, net_size, net_size)
    laplacian = spzeros(Float64, net_size, net_size)
    spectral_history = SpectralHistory([], [], [], [], [], [])
    
    return Network(net_size, positions, adjacency, laplacian, spectral_history)
end


function build_laplacian_sparse!(network::Network)
    adj = network.adjacency
    degrees = vec(sum(adj, dims=2))  # O(n) operation
    
    # Start with -A (reuse sparse structure)
    network.laplacian = -adj

    # Add degrees to diagonal (O(n) operation)
    network.laplacian[diagind(network.laplacian)] .+= degrees
    
    return network.laplacian
end


#The default quantities will be dt = 30min, speed_var = (10m/s)^2 = (0.6 km/min)^2

function simulate!(network::Network, duration_hours::Float64, dt_minutes::Int64 = 30, speed_variance::Float64 = 0.036)
    
    n_steps = Int(duration_hours * 60 / dt_minutes)
    
    for step in 1:n_steps
        # Update positions
        random_walk_step!(network, dt_minutes, speed_variance)
        
        # Rebuild connectivity
        build_adjacency!(network)
        
        # Compute and store spectral properties
        spectral_data = compute_spectral_properties(network)
        update_spectral_history!(network, spectral_data, step * dt_minutes / 60.0)
        
        # Print progress every 10 steps
        if step % 10 == 0
            println("Step $step/$n_steps, Components: $(spectral_data.num_components)")
        end
    end
end


function update_spectral_history!(network::Network, spectral_data, time::Float64)
    push!(network.spectral_history.algebraic_connectivity, spectral_data.algebraic_connectivity)
    push!(network.spectral_history.num_components, spectral_data.num_components)
    
    # Cheeger bounds
    cheeger_lower = 0.5 * spectral_data.spectral_gap
    cheeger_upper = sqrt(2 * spectral_data.max_degree * spectral_data.spectral_gap)
    
    push!(network.spectral_history.cheeger_lower_bound, cheeger_lower)
    push!(network.spectral_history.cheeger_upper_bound, cheeger_upper)
    push!(network.spectral_history.timestamps, time)
end

function spectral_evolution()
    #Create time series of eigenvalue history
end