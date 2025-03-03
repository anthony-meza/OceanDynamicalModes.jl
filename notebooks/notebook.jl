### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c90f0c9e-f5e2-11ef-3763-a79a1560ab8a
begin 
	import Pkg
	Pkg.activate("../")
	Pkg.instantiate()
end

# ╔═╡ 42151ed2-fa30-4f63-b404-bbe19c488a27
begin
	using OceanDynamicalModes
	using LinearAlgebra
	using MeshArrays, PlutoUI
	using Statistics, Plots, StatsPlots

	uniform_N²(z, N²) = N² .* ones(length(z))
end

# ╔═╡ f98e3fe0-db65-47f5-b4cb-2a4f51dddc4e
md"""
Finding internal modes is equivalen to solving the following eigenvalue problem: 

$$\frac{\partial^2 w}{\partial z^2} + c^2 N^2(z) w = 0$$ 

For the case of constant $N^2 = N_0^2$ , the solution os of the form: 

$$w(z) = C_1 \sin(\frac{ N_0 }{c }z) + C_2 \cos(\frac{ N_0 }{c }z)$$ 

Satisifying the top boundary condition ($w(0) = 0$) requires that $C_2$ = 0. The bottom boundary condition ($w(H) = 0$) implies that $\sin(\frac{ N_0 }{c }H) = 0$, which means that $N_0 H / c = n \pi$ for $n = 1, 2, \dots$. 

Linearity and boundary conditions of the problem thus imply that solutions $w_n$ are of the form 

$$\begin{align} w_n(z) &\sim sin(\frac{n \pi}{H} z) \\  c_n &= \frac{N_0 H}{n \pi}\end{align}$$

**Below, we solve for these modes using a uniformly spaces vertical grid and the ECCO V4r4 vertical grid. We compare our results to what is expected from the analytical solution.**
"""

# ╔═╡ 7f6a0978-2071-47dd-8843-cadf02746d80
begin 
    function sol_parameters()
        return PlutoUI.combine() do Child
            # Create an array of actual values (logarithmically spaced)
            N2_values = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]
            
            # Create a selector for actual N² values
            N2_selector = Child("N2", Select(N2_values, default=1e-6))
            
            # For the grid spacing, using a regular slider
            dz_slider = Child("dz", Slider(15:25:500, default=50, show_value=true))
            
            # For the number of modes, using a slider
            nmodes_slider = Child("nmodes", Slider(1:1:6, default=2, show_value=true))
		
            inputs = [
                md""" **N² value:** $(N2_selector) s⁻² """,
                
                md""" **Vertical grid spacing:** $(dz_slider) meters""",
                
                md""" **Number of modes:** $(nmodes_slider) """
            ]
            
            md"""
            #### Uniform Grid and Stratification Parameters
            $(inputs)
            """
        end
    end
    # In another cell, bind the parameters
    @bind strat_params sol_parameters()
end

# ╔═╡ 15f13433-dcd5-4a21-a916-abb84913fa3a
begin 
	#setup ECCO grid and stratification
	pth = MeshArrays.GRID_LLC90
	γ = GridSpec("LatLonCap",pth); Γ = GridLoad(γ;option="full");
	delSq_ECCO = build_delsq_matrix(Γ)
	depths_ECCO = -Γ.RC;
	H = -Γ.RF[end]
	N² = strat_params.N2
	N²_ECCO = diagm(uniform_N²(depths_ECCO, N²))

	#setup uniform grid and stratification for comparison
	dz_uniform = strat_params.dz; 
	depths_uniform = collect(0 + dz_uniform:dz_uniform:H-dz_uniform)
	N²_uniform = diagm(uniform_N²(depths_uniform, N²))
    delSq_uniform = build_delsq_matrix(depths_uniform);
	nothing #done, trick to stop things from printing 
	
end

# ╔═╡ cc490191-fc40-46d9-bedb-1241a04ecd3a
begin 
	nmodes = strat_params.nmodes
	#calculate vertical modes on a uniform grid and ECCO V4r4 vertical gri
	analytical_modes(z, n, H) = @. sin(n * π * z / H)
	
	wmodes_ECCO, ce_ECCO = dynmodes(delSq_ECCO, N²_ECCO; nmodes = nmodes)
	
	wmodes_uni, ce = dynmodes(delSq_uniform, N²_uniform; nmodes = nmodes);

	#calculate analytical modes
	depth_analy = collect(0:1:H)
	wmodes_analy = analytical_modes(depth_analy, collect(1:nmodes)', H);
	nothing #trick to stop things from printing
end

# ╔═╡ 16ae82b0-f6d4-499c-b07e-44cc1d9e1f23
begin 

	titles = ["ECCO V4r4 Grid", "Uniform Grid (Δz = $dz_uniform)", "Analytical Solution"]
	data = [wmodes_ECCO, wmodes_uni, wmodes_analy]
	depths = [depths_ECCO, depths_uniform, depth_analy]
	
	p = plot(layout=(1,3), size=(700, 500), 
	    link=:both,   plot_titlevspan=0.05, 
		plot_title="Vertical Modes and Speeds Comparison (N₀² = $N²)")
	for i in 1:3, j in 1:nmodes
	        plot!(p[i], data[i][:, j], 
	            -depths[i], label="Mode $j", linewidth=2.5)
	    plot!(p[i], title=titles[i])
	    if i == 1
	        plot!(p[i], ylabel="Depth [meters]")
	    end
	end
	
	#speeds for uniform N 
	ce_analytical =  (sqrt(N²) * H) ./ (collect(1:nmodes) .* π)

    # Create descriptive mode labels
    mode_labels = ["Mode $i" for i in 1:nmodes]
    
    # Create grouped bar chart
    p1 = groupedbar(
        mode_labels, 
        [ce_analytical ce_ECCO ce],
        bar_position = :dodge,
        bar_width = 0.7,
        # title = "Wave Speed Comparison by Mode",
        ylabel = "Wave Speed (m/s)",
        legend = :topright,
        label = ["Analytical" "LLC90 Vertical Grid" "Uniformly Spaced Grid"],
        color = [:blue :red :green],
        alpha = 0.8,
        size = (700, 450),
        framestyle = :box,
        fontfamily = "Computer Modern"
    )
    combined_plot = plot(
        p, 
        p1, 
        layout = grid(2, 1, heights=[0.6, 0.4]),
        size = (750, 1000),
    )
end

# ╔═╡ bbc3256c-9ffb-4aad-b4b2-6ad282ee90a2


# ╔═╡ Cell order:
# ╠═c90f0c9e-f5e2-11ef-3763-a79a1560ab8a
# ╠═42151ed2-fa30-4f63-b404-bbe19c488a27
# ╠═15f13433-dcd5-4a21-a916-abb84913fa3a
# ╟─f98e3fe0-db65-47f5-b4cb-2a4f51dddc4e
# ╟─7f6a0978-2071-47dd-8843-cadf02746d80
# ╠═cc490191-fc40-46d9-bedb-1241a04ecd3a
# ╠═16ae82b0-f6d4-499c-b07e-44cc1d9e1f23
# ╠═bbc3256c-9ffb-4aad-b4b2-6ad282ee90a2
