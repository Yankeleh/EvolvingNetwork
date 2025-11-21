#=
Module for computing 
=#


using LinearAlgebra, Arpack, SparseArrays

function compute_spectral_properties(network::Network)
    n::Int64  = network.net_size
    
    # Compute smallest eigenvalues + eigenvectors
    # Request k+1 eigenvalues (0 will be smallest)
    k = min(50, n-2)  # Get 50 smallest eigenvalues
    
    λ, v = eigs(Symmetric(network.laplacian), nev=k, which=:SR, ritzvec=true)
    
    # Get largest eigenvalues (spectral radius) for estimating Cheeger constant
    λ_max, _ = eigs(Symmetric(network.laplacian), nev=2, which=:LR)
    spectral_radius = real(λ_max[1])
    spectral_gap = real(λ_max[1]-λ_max[2])

    # Count zero eigenvalues (connected components)
    zero_tol = 1e-10
    num_components = sum(abs.(λ) .< zero_tol)
    
    # Find algebraic connectivity (smallest non-zero eigenvalue)
    nonzero_idx = findfirst(abs.(λ) .>= zero_tol)
    algebraic_connectivity = nonzero_idx !== nothing ? λ[nonzero_idx] : 0.0
    fiedler_vector = nonzero_idx !== nothing ? v[:, nonzero_idx] : zeros(n)
    
    return (
        eigenvalues = λ,
        spectral_radius = spectral_radius,
        spectral_gap = spectral_gap,
        num_components = num_components,
        algebraic_connectivity = algebraic_connectivity,
        fiedler_vector = fiedler_vector
    )
end