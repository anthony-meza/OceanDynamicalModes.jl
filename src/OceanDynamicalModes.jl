module OceanDynamicalModes

export Grid,
       DynamicalModes

export generate_ECCO_vertical_grid,
       build_delsq_matrix,
       dynmodes,
       compute_dynamical_modes,
       cosine_similarity,
       mode_similarity

include("GriddedModes.jl")

end