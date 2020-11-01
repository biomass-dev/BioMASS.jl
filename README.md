# The BioMASS module for Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://himoto.github.io/BioMASS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://himoto.github.io/BioMASS.jl/dev)
[![Build Status](https://travis-ci.com/himoto/BioMASS.jl.svg?branch=master)](https://travis-ci.com/himoto/BioMASS.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This module provides a Julia interface to the [BioMASS](https://github.com/okadalabipr/biomass) parameter estimation.

## Usage
```julia
const MODEL_PATH = "./fos_model"
include(MODEL_PATH * "/exec_model.jl")
using .ExecModel

using BioMASS

# Estimate unknown model parameters against experimental observations.
optimize(MODEL_PATH, 1, max_generation=20000, allowable_error=0.5)

# Convert optimized parameters into BioMASS format.
param2biomass("./fos_model/")
```

## Installation
```julia
pkg> add https://github.com/himoto/BioMASS.jl  # Press ']' to enter the Pkg REPL mode.
```

## License
[MIT](LICENSE)