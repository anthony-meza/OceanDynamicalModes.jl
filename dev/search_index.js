var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = OceanDynamicalModes","category":"page"},{"location":"#OceanDynamicalModes","page":"Home","title":"OceanDynamicalModes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for OceanDynamicalModes.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [OceanDynamicalModes]","category":"page"},{"location":"#OceanDynamicalModes.build_delsq_matrix-Tuple{Any, Real}","page":"Home","title":"OceanDynamicalModes.build_delsq_matrix","text":"build_delsq_matrix(depth, dz)\n\nConstruct a discretized Laplacian matrix for a given depth and spacing.\n\nThis function generates a tridiagonal matrix representing the Laplacian operator in one dimension for a specific depth and spacing.\n\nArguments\n\ndepth: The depth of the domain.\ndz: The spacing between grid points.\n\nReturns\n\ndelsq: A tridiagonal matrix representing the second derivative operator.\nnz: The number of grid points.\ndz: The grid spacing.\n\n\n\n\n\n","category":"method"},{"location":"#OceanDynamicalModes.clean_up_modes-Tuple{Any, Any, Any}","page":"Home","title":"OceanDynamicalModes.clean_up_modes","text":"clean_up_modes(eigenvalues, wmodes, nmodes)\n\nClean and process eigenvalues and associated modes for further analysis.\n\nThis function takes eigenvalues and their associated modes and      returns a cleaned subset of eigenvalues and modes.\n\nArguments\n\neigenvalues: A 1D array of eigenvalues.\nwmodes: A 2D array where each column represents a mode.\nnmodes: The number of modes to keep.\n\nReturns\n\neigenvalues: A 1D array of cleaned and sorted eigenvalues.\nwmodes: A 2D array where columns are vertical modes \n\n\n\n\n\n","category":"method"},{"location":"#OceanDynamicalModes.dynmodes","page":"Home","title":"OceanDynamicalModes.dynmodes","text":"dynmodes(Nsq, depth, dz, nmodes, rho0=1028)\n\nCalculate dynamic ocean modes \n\nThis function computes dynamic ocean modes, including vertical velocity modes (wmodes), horizontal velocity modes (pmodes), and modal phase speeds (ce) based on a given buoyancy frequency profile (Nsq), depth, and uniform grid spacing (dz).\n\nArguments\n\nNsq: A 1D array representing the squared buoyancy frequency profile.\ndepth: The depth of the ocean.\ndz: The vertical grid spacing.\nnmodes: The number of modes to compute.\nrho0: The reference density of the ocean (default: 1028 kg/m³).\n\nReturns\n\nwmodes: Vertical velocity modes.\npmodes: Horizontal velocity modes.\nce: Modal phase speeds.\n\n\n\n\n\n","category":"function"}]
}
