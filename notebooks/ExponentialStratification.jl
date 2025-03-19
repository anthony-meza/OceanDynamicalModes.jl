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

Below, we solve for the eigenmodes that solve this equation using a uniformly spaces vertical grid and the ECCO V4r4 vertical grid. We use an exponential stratification profile $N(z) = N_0 \exp(z/b)$.
"""

# ╔═╡ 292fb8ed-cbba-4d12-8468-7850848380f1
begin 
    function sol_parameters()
        return PlutoUI.combine() do Child
            # Create an array of actual values (logarithmically spaced)
            N0_values = [1.e-3, 3e-3, 7.3e-3, 8e-3, 9e-3]

            # Create a selector for actual N² values
            N0_selector = Child("N0", Select(N0_values, default=7.3e-3))
            
            # For the grid spacing, using a regular slider
            dz_slider = Child("dz", Slider(15.:25.:500., default=15., show_value=true))
            
            # For the number of modes, using a slider
            nmodes_slider = Child("nmodes", Slider(1:1:6, default=2, show_value=true))
			
			b_slider =  Child("b", Slider(350.:25.:2500., default=1000., show_value=true))

            inputs = [
                md""" **N₀ value:** $(N0_selector) s⁻¹ """,
				md""" **b value:** $(b_slider) meters """,

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
	N₀ = strat_params.N0
	b = strat_params.b
	
	nmodes = strat_params.nmodes
	dz_uniform = strat_params.dz;

	exp_strat(z) = N₀ .* exp.(z / b)
end

# ╔═╡ 7be9c473-12d2-468d-a496-fbf56cc8abee
begin 
	grid_ECCO = generate_ECCO_vertical_grid(z -> exp_strat(z)^2)
	H = grid_ECCO.H

	#setup uniform grid and stratification for comparison
	depths_uniform = -collect(dz_uniform:dz_uniform:H-dz_uniform)
	N_uniform = exp_strat.(depths_uniform)
	N²_uniform = N_uniform.^2
	grid_uniform = Grid(N²_uniform, depths_uniform, dz_uniform, H)
	
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
	
	# Create a consistent color palette for modes
	mode_colors = palette(:tab10, nmodes)
	
	# Create the main plot with 4 subplots in a row
	p = plot(
	    layout = (1, 4),
	    size = (1200, 500),
	    link = :y,         # Link only y-axis for better depth alignment
	    margin = 1Plots.mm,
	    plot_title = "Vertical Modes and Stratification Comparison",
	    plot_titlefontsize = 14,
	    bottom_margin = 5Plots.mm
	)
	
	plot_params = Dict(
    :linewidth => 2.5, :guidefontsize => 10,
    :titlefontsize => 11, :xrotation => 45,
	:legend => :bottomright)

	
	# Panel 1: Stratification Profile
	plot!(
	    p[1], 
	    sqrt.(grid_uniform.N²), 
	    uniform_modes.grid.z, 
		label = nothing,
	    title = "Stratification Profile - N(z)",
	    xlabel = "N [1/s]",
	    ylabel = "Depth [m]",
	    color = :black; plot_params...
	)
	
	# Panel 2: Buoyancy Frequency
	plot!(
	    p[2], 
	    grid_uniform.N², 
	    uniform_modes.grid.z, 
	    label = nothing,
	    title = "Buoyancy Frequency - N²(z)",
	    xlabel = "N² [1/s²]",
	    color = :black; plot_params...

	)
	plot_params[:xlims] = (-1.3, 1.3)
	for j in 1:nmodes
		# Panel 3: ECCO Modes
	    plot!(p[3], ECCO_modes.w[:, j], ECCO_modes.grid.z,
	        label = "Mode $j", 
	        color = mode_colors[j],
			title = "ECCO V4r4 Grid", xlabel = "Amplitude";
			plot_params...
	    )
		# Panel 4: Uniform Modes
		plot!(p[4], uniform_modes.w[:, j], uniform_modes.grid.z,
	        label = "Mode $j", color = mode_colors[j],
			title = "Uniform Grid (Δz = $dz_uniform)",
			xlabel = "Amplitude";
			plot_params...
	    )
		
	end

	# Create wave speed comparison bar chart
	mode_labels = ["Mode $i" for i in 1:nmodes]
	
	p_bar = groupedbar(
	    mode_labels,
	    [ECCO_modes.c uniform_modes.c],
	    bar_position = :dodge,
	    bar_width = 0.7,
	    ylabel = "Wave Speed (m/s)",
	    legend = :topright,
	    label = ["ECCO V4r4 Grid" "Uniform Grid"],
	    color = [mode_colors[1] mode_colors[2]],
	    alpha = 0.8,
	    xtickfontsize = 10,
	    ytickfontsize = 10,
	    guidefontsize = 11,
	    legendfontsize = 9,
	    title = "Wave Speed Comparison",
	    titlefontsize = 12,
	    framestyle = :box, 
	    widen = false,       # Don't widen the plot automatically
	    right_margin = 40Plots.mm,
	    left_margin = 40Plots.mm

	)
	
	# Combine plots with proper sizing
	combined_plot = plot(
	    p, p_bar,
	    layout = grid(2, 1, heights=[0.65, 0.35]),
	    size = (1000, 800),
	    dpi = 300
	)
	
	# Display the plot
	combined_plot
end

# ╔═╡ 4ec396aa-eec4-4504-a8aa-bbfc2f674eb4


# ╔═╡ Cell order:
# ╠═1433b08e-04d6-11f0-09e0-3b3682723ce9
# ╠═b3a26491-f02d-4d01-a2e0-b837cc89fdef
# ╟─522b5767-c6dc-42d0-b6d3-f6f0bcf517b7
# ╟─292fb8ed-cbba-4d12-8468-7850848380f1
# ╟─d39d0089-062d-47b5-bf46-6aa483494f2e
# ╠═7be9c473-12d2-468d-a496-fbf56cc8abee
# ╠═b2d7ae53-99e7-4010-b5bb-ec11e5d04271
# ╠═4ec396aa-eec4-4504-a8aa-bbfc2f674eb4
