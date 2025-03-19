var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = OceanDynamicalModes","category":"page"},{"location":"#OceanDynamicalModes","page":"Home","title":"OceanDynamicalModes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for OceanDynamicalModes.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [OceanDynamicalModes]","category":"page"},{"location":"#OceanDynamicalModes.DynamicalModes","page":"Home","title":"OceanDynamicalModes.DynamicalModes","text":"struct DynamicalModes{T <: Real}\n\nStructure containing dynamical modes information for a stratified fluid.\n\nFields\n\nN²::Union{T,AbstractArray{T,1}}: Buoyancy frequency squared (Brunt-Väisälä frequency squared)\nz::AbstractArray{T,1}: Depth at cell centers\nw::AbstractArray{T,2}: Eigenvectors (vertical modes) for vertical velocity\nv::AbstractArray{T,2}: Eigenvectors for horizontal velocity (if computed)\nc::AbstractArray{T,1}: Modal phase speeds\n\n\n\n\n\n","category":"type"},{"location":"#OceanDynamicalModes.Grid","page":"Home","title":"OceanDynamicalModes.Grid","text":"struct Grid{T <: Real}\n\nGrid structure with vertical discretization information.\n\nFields\n\nΔz::AbstractArray{T,1}: distance between \nz::AbstractArray{T,1}: Depth at cell centers\nnz::Int: Number of vertical grid points\n\n\n\n\n\n","category":"type"},{"location":"#OceanDynamicalModes.build_delsq_matrix-Tuple{Grid}","page":"Home","title":"OceanDynamicalModes.build_delsq_matrix","text":"build_delsq_matrix(g::Grid)\n\nBuild a second-derivative (Laplacian) matrix for the vertical grid. Returns a tridiagonal matrix representing the vertical second derivative operator.\n\n\n\n\n\n","category":"method"},{"location":"#OceanDynamicalModes.compute_dynamical_modes-Union{Tuple{Grid{T}}, Tuple{T}} where T<:AbstractFloat","page":"Home","title":"OceanDynamicalModes.compute_dynamical_modes","text":"compute_dynamical_modes(g::Grid{T}, N²::Union{T,Vector{T}}; nmodes=nothing) where T<:AbstractFloat\n\nCompute and return a DynamicalModes object containing the modes and their properties.\n\nArguments\n\ng::Grid{T}: Grid object containing the vertical discretization\nN²::Union{T,Vector{T}}: Buoyancy frequency squared (Brunt-Väisälä frequency squared)\n\nReturns\n\nmodes::DynamicalModes{T}: Object containing the computed modes and their properties\n\n\n\n\n\n","category":"method"},{"location":"#OceanDynamicalModes.cosine_similarity-Tuple{Any, Any}","page":"Home","title":"OceanDynamicalModes.cosine_similarity","text":"cosine_similarity(v1, v2)\n\nCompute the cosine similarity between two vectors, handling NaN values. Returns a value between 0 and 1, where 1 indicates identical vectors.\n\n\n\n\n\n","category":"method"},{"location":"#OceanDynamicalModes.dynmodes-Union{Tuple{Grid{T}}, Tuple{T}} where T<:AbstractFloat","page":"Home","title":"OceanDynamicalModes.dynmodes","text":"dynmodes(g::Grid{T}, N²::Union{T,Vector{T}}; nmodes=nothing) where T<:AbstractFloat\n\nCompute the dynamical modes (eigenmodes) of a stratified fluid system.\n\nArguments\n\ng::Grid{T}: Grid object containing the vertical discretization\nN²::Union{T,Vector{T}}: Buoyancy frequency squared (Brunt-Väisälä frequency squared)\nCan be a scalar value (constant stratification)\nCan be a vector of length g.nz (depth-dependent stratification)\n\nReturns\n\nwmodes: Matrix of eigenvectors (vertical modes), where each column is a mode\nIf nmodes is specified, returns only the first nmodes modes\nce: Vector of modal phase speeds, computed as 1/√λ where λ are the eigenvalues\nIf nmodes is specified, returns only the first nmodes modal speeds\n\nDetails\n\nThis function solves the generalized eigenvalue problem:\n\ndelsq * w = -λ * N² * w\n\nwhere delsq is the discretized Laplacian operator, N² is the buoyancy frequency squared, and λ are the eigenvalues.\n\n\n\n\n\n","category":"method"},{"location":"#OceanDynamicalModes.generate_ECCO_vertical_grid-Union{Tuple{Union{AbstractVector{T}, T}}, Tuple{T}} where T<:AbstractFloat","page":"Home","title":"OceanDynamicalModes.generate_ECCO_vertical_grid","text":"generate_ECCO_vertical_grid()\n\nGenerate a Grid object using the ECCO vertical grid information. Returns a Grid object with cell thicknesses, depths, and number of levels.\n\n\n\n\n\n","category":"method"},{"location":"#OceanDynamicalModes.mode_similarity-Tuple{DynamicalModes, DynamicalModes}","page":"Home","title":"OceanDynamicalModes.mode_similarity","text":"mode_similarity(modes1::DynamicalModes, modes2::DynamicalModes; \n                nmodes=nothing, interpolate=false, target_grid=1)\n\nCompare eigenvectors using cosine similarity with optional interpolation.\n\nArguments\n\nmodes1, modes2: DynamicalModes objects containing eigenvectors and grid info\nnmodes: Number of modes to compare (default: all)\ninterpolate: Whether to interpolate modes to a common grid\ntarget_grid: Which grid to interpolate onto (1 or 2)\n\nReturns\n\nsimilarity: Matrix of cosine similarity values (values closer to 1 indicate higher similarity)\n\n\n\n\n\n","category":"method"}]
}
