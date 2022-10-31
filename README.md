# The BioMASS module for Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://biomass-dev.github.io/BioMASS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://biomass-dev.github.io/BioMASS.jl/dev)
[![Actions Status](https://github.com/biomass-dev/BioMASS.jl/workflows/CI/badge.svg)](https://github.com/biomass-dev/BioMASS.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Cancers Paper](https://img.shields.io/badge/DOI-10.3390%2Fcancers12102878-blue)](https://doi.org/10.3390/cancers12102878)

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

### Python package requirements:

- numpy - https://numpy.org
- scipy - https://scipy.org
- matplotlib - https://matplotlib.org

## Example

### Model development

This example shows you how to build a simple Michaelis-Menten two-step enzyme catalysis model.

> E + S ⇄ ES → E + P

[`pasmopy.Text2Model`](https://pasmopy.readthedocs.io/en/latest/model_development.html) allows you to build a BioMASS model from text. You simply describe biochemical reactions and the molecular mechanisms extracted from text are converted into an executable model.

Prepare a text file describing the biochemical reactions (e.g., `michaelis_menten.txt`)

```
E + S ⇄ ES | kf=0.003, kr=0.001 | E=100, S=50
ES → E + P | kf=0.002

@obs Substrate: u[S]
@obs E_free: u[E]
@obs E_total: u[E] + u[ES]
@obs Product: u[P]
@obs Complex: u[ES]

@sim tspan: [0, 100]
```

Convert the text into an executable model

```shell
$ python  # pasmopy requires Python 3.7+
```

```python
>>> from pasmopy import Text2Model
>>> description = Text2Model("michaelis_menten.txt", lang="julia")
>>> description.convert()  # generate 'michaelis_menten_jl/'
```

Simulate the model using BioMASS.jl

```shell
$ julia
```

```julia
using BioMASS

model = Model("./michaelis_menten_jl");
run_simulation(model)
```

![michaelis_menten](https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/michaelis_menten_sim.png)

### Parameter estimation

```julia
using BioMASS

model = Model("./examples/fos_model");

# Estimate unknown model parameters from experimental observations
scipy_differential_evolution(model, 1)  # requires scipy package

# Save simulation results to figure/ in the model folder
run_simulation(model, viz_type="best", show_all=true)
```

![estimated_parameter_sets](https://raw.githubusercontent.com/biomass-dev/biomass/master/docs/_static/img/estimated_parameter_sets.png)

## References

- Imoto, H., Zhang, S. & Okada, M. A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data—Application to the ErbB Receptor Signaling Pathway. _Cancers_ **12**, 2878 (2020). https://doi.org/10.3390/cancers12102878

- Imoto, H., Yamashiro, S. & Okada, M. A text-based computational framework for patient -specific modeling for classification of cancers. _iScience_ **25**, 103944 (2022). https://doi.org/10.1016/j.isci.2022.103944

## License

[MIT](https://github.com/biomass-dev/BioMASS.jl/blob/master/LICENSE)
