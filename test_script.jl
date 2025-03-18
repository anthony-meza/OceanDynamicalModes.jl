using OceanDynamicalModes
using PyPlot, LinearAlgebra

dz = 0.05
depth = collect(0 + dz:dz:1. - 0.05)
Nsq = ones(length(depth))
nmodes = 3
wmodes, pmodes, ce = dynmodes(Nsq, depth, dz, nmodes)

ones(19, 3) .* ce'
print(ce)

#compare to analytical solution
c = 1 / pi
analytical = ((c, c / 2, c / 3))

delsq, nz, dz = build_delsq_matrix(depth, dz)
# N-squared diagonal matrix
Nsq_mat = diagm(Nsq)
# Solve generalized eigenvalue problem for eigenvalues and vertical velocity modes
eigenvalues, wmodes = eigen(Array(delsq), -Nsq_mat)
eigenvaluesc, wmodesc = clean_up_modes(eigenvalues, wmodes, 19)

print(analytical)

#plot the modes
fig, ax = plt.subplots()
for i = 1:3
    ax.plot(wmodes[:, i], -depth, label = "Mode " * string(i))
end
ax.legend()
fig


#plot the modes
fig, ax = plt.subplots()
for i = 1:3
    ax.plot(pmodes[:, i], -(depth[2:end] .+ depth[1:end-1]) ./ 2, label = "Mode " * string(i))
end
ax.legend()
fig


