# The BioMASS module for Julia

[![Actions Status](https://github.com/himoto/BioMASS.jl/workflows/CI/badge.svg)](https://github.com/himoto/BioMASS.jl/actions)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://himoto.github.io/BioMASS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://himoto.github.io/BioMASS.jl/dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This module provides a Julia interface to the [BioMASS](https://github.com/okadalabipr/biomass) parameter estimation.

## Usage
```julia
using Distributed
addprocs(10); # add worker processes
@everywhere using BioMASS

@everywhere begin
    model = load_model("./fos_model")
    function optimize_parallel(i)
        optimize(model, i, max_generation=20000, allowable_error=0.5)
    end
end

# Estimate unknown model parameters against experimental observations.
pmap(optimize_parallel, 1:10)

# Save simulation results to figure/ in the model folder
visualize(model, viz_type="best", show_all=true)

# Convert optimized parameters into BioMASS format.
param2biomass("./fos_model")
```

## Installation
```julia
pkg> add https://github.com/himoto/BioMASS.jl  # Press ']' to enter the Pkg REPL mode.
```

## References
- Nakakuki, T. *et al.* Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. *Cell* **141**, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054

- Yao, G., Lee, T. J., Mori, S., Nevins, J. R. & You, L. A bistable Rb-E2F switch underlies the restriction point. *Nat. Cell Biol.* **10**, 476–482 (2008). https://doi.org/10.1038/ncb1711

- Barr, A. R., Heldt, F. S., Zhang, T., Bakal, C. & Novák, B. A Dynamical Framework for the All-or-None G1/S Transition. *Cell Syst.* **2**, 27–37 (2016). https://doi.org/10.1016/j.cels.2016.01.001

- Rata, S. *et al.* Two Interlinked Bistable Switches Govern Mitotic Control in Mammalian Cells. *Curr. Biol.* **28**, 3824-3832.e6 (2018). https://doi.org/10.1016/j.cub.2018.09.059