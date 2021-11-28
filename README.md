# The BioMASS module for Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://biomass-dev.github.io/BioMASS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://biomass-dev.github.io/BioMASS.jl/dev)
[![Actions Status](https://github.com/biomass-dev/BioMASS.jl/workflows/CI/badge.svg)](https://github.com/biomass-dev/BioMASS.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This module provides a Julia interface to the [BioMASS](https://github.com/biomass-dev/biomass) parameter estimation.

## Installation

The package is a registered package, and can be installed with `Pkg.add`.

```julia
julia> using Pkg; Pkg.add("BioMASS")
```

or through the `pkg` REPL mode by typing

```
] add BioMASS
```

## Usage

### Parameter estimation

```julia
using BioMASS

model = Model("./examples/fos_model");

# Estimate unknown model parameters against experimental observations.
optimize(model, 1, max_generation=20000, allowable_error=0.5)

# Save simulation results to figure/ in the model folder
run_simulation(model, viz_type="best", show_all=true)
```

### Conversion of optimized parameters into BioMASS format

```julia
param2biomass("./examples/fos_model")
```

## License

[MIT](https://github.com/biomass-dev/BioMASS.jl/blob/master/LICENSE)
