# OceanDynamicalModes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://anthony-meza.github.io/OceanDynamicalModes.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://anthony-meza.github.io/OceanDynamicalModes.jl/dev/)
[![Build Status](https://github.com/anthony-meza/OceanDynamicalModes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anthony-meza/OceanDynamicalModes.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anthony-meza/OceanDynamicalModes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anthony-meza/OceanDynamicalModes.jl)

This code solves the generalized eigenvalue problem of the form 
$$\\frac{\\partial^2}{\\partial z^2} w_m + \\alpha^2 N^2 w_m = 0$$
where $w_m = 0$ at the surface and the bottom


Code is based on this Python Code: 

https://gist.github.com/douglatornell/5479638
