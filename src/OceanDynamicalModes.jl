module OceanDynamicalModes
using LinearAlgebra
export build_delsq_matrix, dynmodes, clean_up_modes
"""
    build_delsq_matrix(depth, dz)

Construct a discretized Laplacian matrix for a given depth and spacing.

This function generates a tridiagonal matrix representing the Laplacian operator
in one dimension for a specific depth and spacing.

# Arguments
- `depth`: The depth of the domain.
- `dz`: The spacing between grid points.

# Returns
- `delsq`: A tridiagonal matrix representing the second derivative operator.
- `nz`: The number of grid points.
- `dz`: The grid spacing.
"""
function build_delsq_matrix(depth, dz::Real)
    zed = -depth
    nz = length(zed)    
    d = -2/(dz^2) .* ones(nz)
    dl = 1/(dz^2) .* ones(nz-1)
    du = 1/(dz^2) .* ones(nz-1)
    delsq = Tridiagonal(dl, d, du)
    return delsq, nz, dz
end

"""
    dynmodes(Nsq, depth, dz, nmodes, rho0=1028)

Calculate dynamic ocean modes 

This function computes dynamic ocean modes, including vertical velocity modes (wmodes), horizontal velocity modes (pmodes),
and modal phase speeds (ce) based on a given buoyancy frequency profile (Nsq), depth, and uniform grid spacing (dz).

# Arguments
- `Nsq`: A 1D array representing the squared buoyancy frequency profile.
- `depth`: The depth of the ocean.
- `dz`: The vertical grid spacing.
- `nmodes`: The number of modes to compute.
- `rho0`: The reference density of the ocean (default: 1028 kg/m³).

# Returns
- `wmodes`: Vertical velocity modes.
- `pmodes`: Horizontal velocity modes.
- `ce`: Modal phase speeds.
"""
function dynmodes(Nsq, depth, dz, nmodes, rho0=1028)
    # Del-squared matrix plus boundary conditions
    delsq, nz, dz = build_delsq_matrix(depth, dz)
    # N-squared diagonal matrix
    Nsq_mat = diagm(Nsq)
    # Solve generalized eigenvalue problem for eigenvalues and vertical velocity modes
    eigenvalues, wmodes = eigen(Array(delsq), -Nsq_mat)
    eigenvalues, wmodes = clean_up_modes(eigenvalues, wmodes, nmodes)
    # Modal speeds
    ce = 1.0 ./ sqrt.(real.(eigenvalues))

    # Horizontal velocity modes (NOT COMPUTED YET)
    pmodes = similar(wmodes)
    # 1st derivative of vertical modes
    # pr = diff(wmodes, dims=2) .* rho0 .* (ce .^ 2) .* dz'
    # Linear interpolation on depth grid
    # pmodes[:, 2:nz] .= (pr[:, 2:nz] + pr[:, 1:nz-1]) / 2
    # pmodes[:, 1] .= pr[:, 1]
    # pmodes[:, nz] .= pr[:, nz-1]
    return wmodes, pmodes, ce
end

"""
    clean_up_modes(eigenvalues, wmodes, nmodes)

Clean and process eigenvalues and associated modes for further analysis.

This function takes eigenvalues and their associated modes, applies several filters,
and sorts them for subsequent analysis. It returns a cleaned subset of eigenvalues and modes.

# Arguments
- `eigenvalues`: A 1D array of eigenvalues.
- `wmodes`: A 2D array where each column represents a mode.
- `nmodes`: The number of modes to keep.

# Returns
- `eigenvalues`: A 1D array of cleaned and sorted eigenvalues.
- `wmodes`: A 2D array of corresponding modes.
"""

function clean_up_modes(eigenvalues, wmodes, nmodes)
    # Transpose modes to be handled as an array of vectors
    wmodes = transpose(wmodes)
    
    # Step 1: Filter out complex-valued eigenvalues
    # Check if the imaginary part of eigenvalues is zero
    mask = imag.(eigenvalues) .== 0
    eigenvalues = eigenvalues[mask]  # Keep only real eigenvalues
    wmodes = wmodes[mask, :]  # Corresponding modes
    
    # Step 2: Filter out small/negative eigenvalues
    mask = eigenvalues .>= 1e-10  # Select eigenvalues >= 1e-10
    eigenvalues = eigenvalues[mask]
    wmodes = wmodes[mask, :]
    
    # Step 3: Sort eigenvalues and modes and truncate to the number of modes requested
    index = sortperm(eigenvalues)  # Get the indices for sorting eigenvalues
    eigenvalues = eigenvalues[index[1:nmodes]]  # Sort and keep top nmodes eigenvalues
    wmodes = wmodes[index[1:nmodes], :]  # Sort modes accordingly, each mode corresponds to a particular row
    
    # Step 4: Return the cleaned-up eigenvalues and modes
    return eigenvalues, wmodes
end

end