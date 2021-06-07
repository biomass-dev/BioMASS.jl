# The BioMASS module for Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://biomass-dev.github.io/BioMASS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://biomass-dev.github.io/BioMASS.jl/dev)
[![Actions Status](https://github.com/biomass-dev/BioMASS.jl/workflows/CI/badge.svg)](https://github.com/biomass-dev/BioMASS.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This module provides a Julia interface to the [BioMASS](https://github.com/okadalabipr/biomass) parameter estimation.

![](docs/src/assets/result.png)

## Features

BioMASS.jl supports:

- parameter estimation of ODE/DDE models
- visualization of simulation results
- bifurcation analysis

## Usage

```julia
using BioMASS

model = load_model("./examples/fos_model");

# Estimate unknown model parameters against experimental observations.
optimize(model, 1, max_generation=20000, allowable_error=0.5)

# Save simulation results to figure/ in the model folder
run_simulation(model, viz_type="best", show_all=true)
```

### Conversion of optimized parameters into BioMASS format

![](docs/src/assets/conversion.png)

```julia
param2biomass("./examples/fos_model")
```

## Installation

The package is a registered package, and can be installed with `Pkg.add`.

```julia
julia> using Pkg; Pkg.add("BioMASS")
```

or through the `pkg` REPL mode by typing

```
] add BioMASS
```

## References

- Nakakuki, T. _et al._ Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. _Cell_ **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

- Inoue, K. _et al._ Oscillation dynamics underlie functional switching of NF-κB for B-cell activation. _npj Syst. Biol. Appl._ **2**, 16024 (2016). https://doi.org/10.1038/npjsba.2016.24

- Yao, G., Lee, T. J., Mori, S., Nevins, J. R. & You, L. A bistable Rb-E2F switch underlies the restriction point. _Nat. Cell Biol._ **10**, 476–482 (2008). https://doi.org/10.1038/ncb1711

- Barr, A. R., Heldt, F. S., Zhang, T., Bakal, C. & Novák, B. A Dynamical Framework for the All-or-None G1/S Transition. _Cell Syst._ **2**, 27–37 (2016). https://doi.org/10.1016/j.cels.2016.01.001

- Rata, S. _et al._ Two Interlinked Bistable Switches Govern Mitotic Control in Mammalian Cells. _Curr. Biol._ **28**, 3824-3832.e6 (2018). https://doi.org/10.1016/j.cub.2018.09.059

## License

[MIT](https://github.com/biomass-dev/BioMASS.jl/blob/master/LICENSE)
