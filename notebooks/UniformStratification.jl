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

# ╔═╡ 1433b08e-04d6-11f0-09e0-3b3682723ce9
begin 
	import Pkg
	Pkg.activate("../")
	Pkg.instantiate()
end

# ╔═╡ b3a26491-f02d-4d01-a2e0-b837cc89fdef
begin
	using OceanDynamicalModes
	using LinearAlgebra
	using MeshArrays, PlutoUI
	using Statistics, Plots, StatsPlots
end

# ╔═╡ 522b5767-c6dc-42d0-b6d3-f6f0bcf517b7
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

# ╔═╡ 292fb8ed-cbba-4d12-8468-7850848380f1
begin 
    function sol_parameters()
        return PlutoUI.combine() do Child
            # Create an array of actual values (logarithmically spaced)
            N2_values = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]
            
            # Create a selector for actual N² values
            N2_selector = Child("N2", Select(N2_values, default=1e-6))
            
            # For the grid spacing, using a regular slider
            dz_slider = Child("dz", Slider(15.:25.:500., default=50., show_value=true))
            
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

# ╔═╡ d39d0089-062d-47b5-bf46-6aa483494f2e
begin 
	N² = strat_params.N2
	nmodes = strat_params.nmodes
	dz_uniform = strat_params.dz;
	
	grid_ECCO = generate_ECCO_vertical_grid(N²)
	H = grid_ECCO.H

	#setup uniform grid and stratification for comparison
	depths_uniform = -collect(0 + dz_uniform:dz_uniform:H-dz_uniform)
	grid_uniform = Grid(N², depths_uniform, dz_uniform, H)
	
	#calculate vertical modes on a uniform grid and ECCO V4r4 vertical gri
	ECCO_modes = compute_dynamical_modes(grid_ECCO; nmodes = nmodes)	
	uniform_modes = compute_dynamical_modes(grid_uniform; nmodes = nmodes)	
	
	#calculate analytical modes
	analytical_modes(z, n, H) = @. sin(n * π * z / H)
	depth_analy = -collect(0:1:H)
	wmodes_analy = analytical_modes(depth_analy, collect(1:nmodes)', H);
	nothing #trick to stop things from printing
end

# ╔═╡ b2d7ae53-99e7-4010-b5bb-ec11e5d04271
begin 

	titles = ["ECCO V4r4 Grid", "Uniform Grid (Δz = $dz_uniform)", "Analytical Solution"]
	data = [ECCO_modes.w, uniform_modes.w, wmodes_analy]
	depths = [ECCO_modes.grid.z, uniform_modes.grid.z, depth_analy]
	
	p = plot(layout=(1,3), size=(700, 500), 
	    link=:both,   plot_titlevspan=0.05, 
		plot_title="Vertical Modes and Speeds Comparison (N₀² = $N²)")
	for i in 1:3, j in 1:nmodes
	        plot!(p[i], data[i][:, j], 
	            depths[i], label="Mode $j", linewidth=2.5)
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
        [ce_analytical ECCO_modes.c uniform_modes.c],
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

# ╔═╡ Cell order:
# ╠═1433b08e-04d6-11f0-09e0-3b3682723ce9
# ╠═b3a26491-f02d-4d01-a2e0-b837cc89fdef
# ╟─522b5767-c6dc-42d0-b6d3-f6f0bcf517b7
# ╟─292fb8ed-cbba-4d12-8468-7850848380f1
# ╠═d39d0089-062d-47b5-bf46-6aa483494f2e
# ╠═b2d7ae53-99e7-4010-b5bb-ec11e5d04271
