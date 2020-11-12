# The BioMASS module for Julia

[![Actions Status](https://github.com/himoto/BioMASS.jl/workflows/CI/badge.svg)](https://github.com/himoto/BioMASS.jl/actions)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://himoto.github.io/BioMASS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://himoto.github.io/BioMASS.jl/dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This module provides a Julia interface to the [BioMASS](https://github.com/okadalabipr/biomass) parameter estimation.

## Usage
```julia
using BioMASS

model = load_model("./fos_model")

# Estimate unknown model parameters against experimental observations.
optimize(model, 1, max_generation=20000, allowable_error=0.5)

# Save simulation results to figure/ in the model folder
visualize(model, viz_type="best", show_all=true)

# Convert optimized parameters into BioMASS format.
param2biomass("./fos_model")
```

## Installation
```julia
pkg> add https://github.com/himoto/BioMASS.jl  # Press ']' to enter the Pkg REPL mode.
```