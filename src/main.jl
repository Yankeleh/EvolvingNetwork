#= 
Main module file for EvolvingNetork
Defines module structure, includes submodules, and exports public API
=#

module EvolvingNetwork

using LinearAlgebra, SparseArrays, Random
using StatsBase, Printf
using Geodesy, Arpack
using Base.Threads

# Set random seed for reproducibility
Random.seed!(314)

# Include submodules
include("core.jl")
include("connectivity.jl")
include("dynamics.jl")
include("spectral.jl")
include("analysis.jl")
include("plotting.jl")

# Public API
export Network, 
       SpectralHistory,
       simulate!,
       analyze_connectivity_probability,
       analyze_average_components,
       plot_connectivity_probability,
       plot_average_components,
       plot_cheeger_bounds
end # module