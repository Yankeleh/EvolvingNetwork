#=
Module for running random walk
=#

using Random


#Will use dt = 30min, speed_var = (10m/s)^2 = (0.6 km/min)^2


function random_walk_step!(network::Network, dt::Float64, speed_variance::Float64)
    for i in 1:network.net_size
        # Current position
        lat, lon, alt = network.positions[i, :]
        
        # Sample displacement in kilometers
        Δx = randn() * sqrt(speed_variance * dt)  # East-West displacement
        Δy = randn() * sqrt(speed_variance * dt)  # North-South displacement
        
        # Convert km displacement to degrees (approximate)
        Δlat = Δy / 111.2  # ~111.2 km per degree latitude
        Δlon = Δx / (111.2 * cosd(lat))  # Adjust for longitude convergence
        
        # Update position with cutoff
        max_displacement = 3 * sqrt(speed_variance * dt)
        if sqrt(Δx^2 + Δy^2) <= max_displacement
            network.positions[i, 1] += Δlat
            network.positions[i, 2] += Δlon
            # Keep altitude constant
        end
    end
end