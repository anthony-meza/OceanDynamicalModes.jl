### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8599b2b6-729c-11ee-0253-f34553e62bfa
import Pkg; Pkg.activate("."); Pkg.instantiate()

# ╔═╡ 3c2a77e6-ffde-40a5-b4a6-dd44a9088954
using Pluto, PlutoUI

# ╔═╡ 22255253-1b42-4a79-a56e-4844dddf9656
using OceanDynamicalModes; using PythonPlot

# ╔═╡ 1f0a2ea8-327f-416e-bab2-bfd62c0aa7ac
vertical_structure = @bind vstructure Select(["uniform","exponential"])

# ╔═╡ d7e3e970-d28f-4cc9-a8f0-061df796e0fd
dz = 0.05; depth = collect(0 + dz:dz:1. - 0.05)
	

# ╔═╡ 38ef889a-a5e9-4b87-b70d-fa1801d1bc9e
vstructure == "uniform" ? Nsq = ones(length(depth)) : println("not implemented")

# ╔═╡ a0bdb328-d2dc-42d3-b79c-7bd8016bbeab
begin

	nmodes = 3
	wmodes, pmodes, ce = dynmodes(Nsq, depth, dz, nmodes)
	ones(19, 3) .* ce'
	print(ce)

	#compare to analytical solution
	c = 1 / pi
	analytical = ((c, c / 2, c / 3))
	print(analytical)
end

# ╔═╡ 5362e727-8708-4797-9409-858d3ee9e757
begin
	#plot the modes
	fig, ax = subplots()
	for i = 1:3
    	ax.plot(wmodes[:, i], -depth, label = "Mode " * string(i))
	end
	ax.legend()
	fig
end

# ╔═╡ d49994ca-e3d7-4820-9698-acd7b232fb94
begin
	let
#plot the modes
fig, ax = subplots()
for i = 1:3
    ax.plot(pmodes[:, i], -(depth[2:end] .+ depth[1:end-1]) ./ 2, label = "Mode " * string(i))
end
ax.legend()
fig
	end
	
end

# ╔═╡ ebfa3bf8-2198-4026-8cd3-3719373c86fa
a_slider = @bind a Slider(0:.1:100,show_value=true,default=10)

# ╔═╡ ebd5180b-a343-4fed-8bc1-6e05cdc2d64d
println(a)

# ╔═╡ Cell order:
# ╠═8599b2b6-729c-11ee-0253-f34553e62bfa
# ╠═3c2a77e6-ffde-40a5-b4a6-dd44a9088954
# ╠═22255253-1b42-4a79-a56e-4844dddf9656
# ╠═1f0a2ea8-327f-416e-bab2-bfd62c0aa7ac
# ╠═d7e3e970-d28f-4cc9-a8f0-061df796e0fd
# ╠═38ef889a-a5e9-4b87-b70d-fa1801d1bc9e
# ╠═a0bdb328-d2dc-42d3-b79c-7bd8016bbeab
# ╠═5362e727-8708-4797-9409-858d3ee9e757
# ╠═d49994ca-e3d7-4820-9698-acd7b232fb94
# ╠═ebfa3bf8-2198-4026-8cd3-3719373c86fa
# ╠═ebd5180b-a343-4fed-8bc1-6e05cdc2d64d
