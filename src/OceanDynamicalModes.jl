module OceanDynamicalModes
using LinearAlgebra, MeshArrays
export build_delsq_matrix, dynmodes

function build_delsq_matrix(depth::Vector)
	dRC = vcat(0, abs.(depth))
	dRC = dRC[2:end] - dRC[1:end-1]
	dRC = abs.(vcat(dRC, dRC[end] / 2))


    d = -2 * inv.(dRC[1:end-1] .* dRC[2:end])
    dl = 2 * inv.(dRC[2:end-1] .* (dRC[2:end-1] .+ dRC[3:end]))
    du = 2 * inv.(dRC[2:end-1] .* (dRC[2:end-1] .+ dRC[1:end-2]))

    delsq = Matrix(Tridiagonal(dl, d, du))
    return delsq
end

"""
Building second derivative operator on the LLC90 vertical grid
"""
function build_delsq_matrix(Γ::NamedTuple)
    # dRF = abs.(Γ.DRF) #includes an upper ghost point
    dRC = abs.(Γ.DRC);
    #assume that points are spaced equally far away 
    dRC = vcat(dRC, dRC[end] / 2)
    
    d = -2 * inv.(dRC[1:end-1] .* dRC[2:end])
    dl = 2 * inv.(dRC[2:end-1] .* (dRC[2:end-1] .+ dRC[3:end]))
    du = 2 * inv.(dRC[2:end-1] .* (dRC[2:end-1] .+ dRC[1:end-2]))

    delSq = Matrix(Tridiagonal(dl, d, du))
    return delSq
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
function dynmodes(delsq, Nsq_mat; nmodes = nothing)
    eigenvalues, wmodes = eigen(Array(delsq), -Nsq_mat)
    eigenvalues = real.(eigenvalues)

    wmodes = wmodes .* sign.(eigenvalues)
    eigenvalues = eigenvalues .* sign.(eigenvalues)

    # Modal speeds
    ce = 1.0 ./ sqrt.(real.(eigenvalues))

    if ~isnothing(nmodes)
        return wmodes[:, 1:nmodes], ce[1:nmodes]
    else
        return wmodes, ce
    end
end

end
