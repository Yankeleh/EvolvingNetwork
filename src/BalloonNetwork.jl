#= 
Module to implement balloon Network structure. This is what we will evolve.
=#

module BalloonNetwork

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


#Balloon network structure

mutable struct Network
    n_balloons::Int
    positions::Matrix{Float64} #n_balloons x 3 (lat_deg, lon_deg, alt_km)
    
    adjacency::SparseMatrixCSC{Float64, Int}
    laplacian::SparseMatrixCSC{Float64, Int}
    
    spectral_history::SpectralHistory
end

include("dynamics.jl") 
include("connectivity.jl")
include("spectral.jl")
include("analysis.jl")    


#Constructor for the network


function Network(n_balloons::Int)
    # Initialize positions using proper sphere sampling

    positions = Matrix{Float64}(undef, n_balloons, 3)

    for i in 1:n_balloons
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
    adjacency = spzeros(Float64, n_balloons, n_balloons)
    laplacian = spzeros(Float64, n_balloons, n_balloons)
    spectral_history = SpectralHistory([], [], [], [], [], [])
    
    return Network(n_balloons, positions, adjacency, laplacian, spectral_history)
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


function simulate!()
    #Run dynamics.jl
    #Update Network struct
end


function analyze_spectral_evolution()
    #Create time series of eigenvalue history
end




export Network, simulate!, analyze_spectral_evolution



end