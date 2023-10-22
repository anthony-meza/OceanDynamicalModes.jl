module OceanDynamicalModes
using LinearAlgebra

function build_delsq_matrix(depth)
    zed = -depth
    # Depth steps
    nz = length(zed)
    dz = zed[1:nz-1] .- zed[2:nz]
    # Depth step midpoints
    z_mid = zed[1:nz-1] .- 0.5 .* dz
    # Depth step midpoint steps
    dz_mid = similar(zed)
    dz_mid[2:nz-1] .= z_mid[1:nz-2] .- z_mid[2:nz-1]
    dz_mid[1] = dz_mid[2]
    dz_mid[nz] = dz_mid[nz-1]
    # del-squared matrix plus boundary conditions
    delsq = zeros(Float64, nz, nz)
    diag_index = 2:nz-1
    delsq[diag_index, diag_index] .= 1.0 ./ (dz[1:nz-2] .* dz_mid[2:nz-1]) .+ 1.0 ./ (dz[2:nz-1] .* dz_mid[2:nz-1])
    diag_index_p1 = 3:nz
    delsq[diag_index, diag_index_p1] .= -1.0 ./ (dz[1:nz-2] .* dz_mid[2:nz-1])
    diag_index_m1 = 1:nz-2
    delsq[diag_index, diag_index_m1] .= -1.0 ./ (dz[2:nz-1] .* dz_mid[2:nz-1])
    # Boundary conditions
    delsq[1, 1] = delsq[nz, 1] = -1.0
    return delsq, nz, dz
end

function clean_up_modes(eigenvalues, wmodes, nmodes)
    # Transpose modes to be handled as an array of vectors
    wmodes = transpose(wmodes)
    # Filter out complex-values and small/negative eigenvalues and corresponding modes
    mask = imag.(eigenvalues) .== 0
    eigenvalues = eigenvalues[mask]
    wmodes = wmodes[mask]
    mask = eigenvalues .>= 1e-10
    eigenvalues = eigenvalues[mask]
    wmodes = wmodes[mask]
    # Sort eigenvalues and modes and truncate to the number of modes requested
    index = sortperm(eigenvalues)
    eigenvalues = eigenvalues[index[1:nmodes]]
    wmodes = wmodes[index[1:nmodes]]
    return eigenvalues, wmodes
end


function dynmodes(Nsq, depth, nmodes, rho0=1028)
    # Del-squared matrix plus boundary conditions
    delsq, nz, dz = build_delsq_matrix(depth)
    # N-squared diagonal matrix
    Nsq_mat = Matrix(1.0I, nz, nz) * Nsq
    # Solve generalized eigenvalue problem for eigenvalues and vertical velocity modes
    eigenvalues, wmodes = eigen(delsq, Nsq_mat)
    eigenvalues, wmodes = clean_up_modes(eigenvalues, wmodes, nmodes)
    # Modal speeds
    ce = 1.0 ./ sqrt.(real.(eigenvalues))
    # Horizontal velocity modes
    pmodes = similar(wmodes)
    # 1st derivative of vertical modes
    pr = diff(wmodes, dims=2) .* rho0 .* (ce .^ 2) .* dz'
    # Linear interpolation on depth grid
    pmodes[:, 2:nz] .= (pr[:, 2:nz] + pr[:, 1:nz-1]) / 2
    pmodes[:, 1] .= pr[:, 1]
    pmodes[:, nz] .= pr[:, nz-1]
    return wmodes, pmodes, ce
end

end
