using LinearAlgebra
using Interpolations
using MeshArrays



"""
    struct Grid{T <: Real}

Grid structure with vertical discretization information.
# Fields
- `Δz::AbstractArray{T,1}`: distance between 
- `z::AbstractArray{T,1}`: Depth at cell centers
- `nz::Int`: Number of vertical grid points
"""
struct Grid{T <: Real}
    N²::Union{T,AbstractArray{T,1}}
    z::AbstractArray{T,1}
    Δz::Union{T,AbstractArray{T,1}}
    nz::Int
    H::T

    # Inner constructor that automatically sets nz
    function Grid(N²::Union{T,AbstractArray{T,1}}, z::AbstractArray{T,1}, 
        Δz::Union{T,AbstractArray{T,1}}, H::T) where T <: Real
        nz = length(z)

        # Convert Δz to a vector if it's a scalar
        Δz_vec = if isa(Δz, Number)
            fill(Δz, nz + 1)  # Create a vector filled with the scalar value
        else
            1 * Δz  # Already a vector, no change needed
        end
        
        new{T}(N², z, Δz_vec, nz, H)
    end
end

"""
    struct DynamicalModes{T <: Real}

Structure containing dynamical modes information for a stratified fluid.
# Fields
- `N²::Union{T,AbstractArray{T,1}}`: Buoyancy frequency squared (Brunt-Väisälä frequency squared)
- `z::AbstractArray{T,1}`: Depth at cell centers
- `w::AbstractArray{T,2}`: Eigenvectors (vertical modes) for vertical velocity
- `v::AbstractArray{T,2}`: Eigenvectors for horizontal velocity (if computed)
- `c::AbstractArray{T,1}`: Modal phase speeds
"""
struct DynamicalModes{T <: Real}
    grid::Grid
    w::AbstractArray{T,2}
    v::AbstractArray{T,2}
    c::AbstractArray{T,1}
end

"""
    generate_ECCO_vertical_grid()

Generate a Grid object using the ECCO vertical grid information.
Returns a Grid object with cell thicknesses, depths, and number of levels.
"""
function generate_ECCO_vertical_grid(N²::Union{T,AbstractArray{T,1}})  where T<:AbstractFloat
    pth = MeshArrays.GRID_LLC90
    γ = GridSpec("LatLonCap", pth); Γ = GridLoad(γ; option="full");
    Δz = abs.(Γ.DRC);
    #assume that points are spaced equally far away 
    Δz = vcat(Δz, Δz[end] / 2)
    depths = -abs.(Γ.RC);
    H = -Γ.RF[end]
    return Grid(N², depths, Δz, H)
end

function generate_ECCO_vertical_grid(N²_func::Function)
    pth = MeshArrays.GRID_LLC90
    γ = GridSpec("LatLonCap", pth); Γ = GridLoad(γ; option="full");
    Δz = abs.(Γ.DRC);
    Δz = vcat(Δz, Δz[end] / 2)
    depths = -abs.(Γ.RC);
    
    # Apply the function to generate N² values at each depth
    N² = N²_func.(depths)
    
    H = -Γ.RF[end]
    return Grid(N², depths, Δz, H)
end

"""
    build_delsq_matrix(g::Grid)

Build a second-derivative (Laplacian) matrix for the vertical grid.
Returns a tridiagonal matrix representing the vertical second derivative operator.
"""
function build_delsq_matrix(g::Grid)
	dRC = g.Δz

    d = -2 * inv.(dRC[1:end-1] .* dRC[2:end])
    dl = 2 * inv.(dRC[2:end-1] .* (dRC[2:end-1] .+ dRC[3:end]))
    du = 2 * inv.(dRC[2:end-1] .* (dRC[2:end-1] .+ dRC[1:end-2]))

    delsq = Matrix(Tridiagonal(dl, d, du))
    return delsq
end

"""
    dynmodes(g::Grid{T}, N²::Union{T,Vector{T}}; nmodes=nothing) where T<:AbstractFloat

Compute the dynamical modes (eigenmodes) of a stratified fluid system.

# Arguments
- `g::Grid{T}`: Grid object containing the vertical discretization
- `N²::Union{T,Vector{T}}`: Buoyancy frequency squared (Brunt-Väisälä frequency squared)
  - Can be a scalar value (constant stratification)
  - Can be a vector of length `g.nz` (depth-dependent stratification)

# Returns
- `wmodes`: Matrix of eigenvectors (vertical modes), where each column is a mode
  - If `nmodes` is specified, returns only the first `nmodes` modes
- `ce`: Vector of modal phase speeds, computed as `1/√λ` where `λ` are the eigenvalues
  - If `nmodes` is specified, returns only the first `nmodes` modal speeds

# Details
This function solves the generalized eigenvalue problem:
```
delsq * w = -λ * N² * w
```
where `delsq` is the discretized Laplacian operator, `N²` is the buoyancy frequency squared,
and `λ` are the eigenvalues.
"""
function dynmodes(g::Grid{T}; 
                  nmodes = nothing) where T<:AbstractFloat
    # Build the Laplacian matrix
    delsq = build_delsq_matrix(g)
    
    # Create buoyancy frequency matrix
    Nsq_vec = zeros(g.nz)
    Nsq_vec .= g.N²
    Nsq_mat = Matrix(Diagonal(Nsq_vec))
    
    # Solve the eigenvalue problem
    eigenvalues, wmodes = eigen(Array(delsq), -Nsq_mat)
    eigenvalues = real.(eigenvalues)
    
    # Ensure eigenvalues are always positive
    sign_correction = sign.(eigenvalues)
    wmodes .*= sign_correction
    eigenvalues .*= sign_correction
    
    # Modal speeds
    ce = 1.0 ./ sqrt.(max.(eigenvalues, eps(T)))
    
    # Return the requested number of modes
    return isnothing(nmodes) ? (wmodes, ce) : (wmodes[:, 1:nmodes], ce[1:nmodes])
end

"""
    compute_dynamical_modes(g::Grid{T}, N²::Union{T,Vector{T}}; nmodes=nothing) where T<:AbstractFloat

Compute and return a DynamicalModes object containing the modes and their properties.

# Arguments
- `g::Grid{T}`: Grid object containing the vertical discretization
- `N²::Union{T,Vector{T}}`: Buoyancy frequency squared (Brunt-Väisälä frequency squared)

# Returns
- `modes::DynamicalModes{T}`: Object containing the computed modes and their properties
"""
function compute_dynamical_modes(g::Grid{T}; 
                               nmodes = nothing, 
                               calculate_horizontal_modes = false) where T<:AbstractFloat
    # Compute the dynamical modes
    wmodes, ce = dynmodes(g, nmodes=nmodes)
    
    # Compute horizontal velocity modes (v) from vertical velocity modes (w)
    # This is a simplified approach; exact implementation depends on your equations
    vmodes = similar(wmodes)

    if calculate_horizontal_modes
        # for i in 1:size(wmodes, 2)
        #     # Calculate horizontal velocity modes from vertical velocity modes
        #     # This is just an example; modify as needed for your specific problem
        #     vmodes[:, i] = -diff(vcat(0, wmodes[:, i])) ./ g.Δz
        # end
        nothing
    else
        vmodes .= NaN
    end
    
    return DynamicalModes(g, wmodes, vmodes, ce)
end

"""
    cosine_similarity(v1, v2)

Compute the cosine similarity between two vectors, handling NaN values.
Returns a value between 0 and 1, where 1 indicates identical vectors.
"""
function cosine_similarity(v1, v2)
    # Find indices where both vectors have valid values
    valid_idx = .!isnan.(v1) .& .!isnan.(v2)
    
    # Return 0 if no valid points for comparison
    if !any(valid_idx)
        return 0.0
    end
    
    # Use only valid data points for similarity calculation
    v1_valid = v1[valid_idx]
    v2_valid = v2[valid_idx]
    
    # Handle zero-norm vectors
    n1 = norm(v1_valid)
    n2 = norm(v2_valid)
    
    if n1 < eps(Float64) || n2 < eps(Float64)
        return 0.0
    end
    
    # Account for sign flips with abs()
    return abs(dot(v1_valid, v2_valid) / (n1 * n2))
end

"""
    mode_similarity(modes1::DynamicalModes, modes2::DynamicalModes; 
                    nmodes=nothing, interpolate=false, target_grid=1)

Compare eigenvectors using cosine similarity with optional interpolation.

# Arguments
- `modes1`, `modes2`: DynamicalModes objects containing eigenvectors and grid info
- `nmodes`: Number of modes to compare (default: all)
- `interpolate`: Whether to interpolate modes to a common grid
- `target_grid`: Which grid to interpolate onto (1 or 2)

# Returns
- `similarity`: Matrix of cosine similarity values (values closer to 1 indicate higher similarity)
"""
function mode_similarity(modes1::DynamicalModes, modes2::DynamicalModes; 
                         nmodes=nothing, interpolate=false, target_grid=1)
    # Get modes to compare
    n1 = isnothing(nmodes) ? size(modes1.w, 2) : min(nmodes, size(modes1.w, 2))
    n2 = isnothing(nmodes) ? size(modes2.w, 2) : min(nmodes, size(modes2.w, 2))
    
    # Set up modes for comparison
    m1 = modes1.w[:, 1:n1]
    m2 = modes2.w[:, 1:n2]
    
    # Interpolate if requested
    if interpolate
        
        if target_grid == 1
            # Interpolate modes2 onto grid1
            m2_interp = similar(m1)
            for j in 1:n2
                itp = LinearInterpolation(modes2.z, m2[:, j], extrapolation_bc=NaN)
                m2_interp[:, j] = itp.(modes1.z)
            end
            m2 = m2_interp
        else
            # Interpolate modes1 onto grid2
            m1_interp = similar(m2)
            for i in 1:n1
                itp = LinearInterpolation(modes1.z, m1[:, i], extrapolation_bc=NaN)
                m1_interp[:, i] = itp.(modes2.z)
            end
            m1 = m1_interp
        end
    end
    
    # Calculate cosine similarity
    similarity = zeros(eltype(m1), n1, n2)
    for i in 1:n1
        for j in 1:n2
            # Get vectors to compare
            v1 = m1[:,i]
            v2 = m2[:,j]
            
            # Use separate cosine similarity function
            similarity[i, j] = cosine_similarity(v1, v2)
        end
    end
    
    return similarity
end