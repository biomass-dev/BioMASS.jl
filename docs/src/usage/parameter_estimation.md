# Parameter Estimation

**load_model**(```path_to_model```::String)
- Loads a BioMASS model. The model must include the following files:

|Name|Content|
|:--|:--|
|`name2idx/`|Names of model parameters and species|
|`set_model.jl`|Differential equation, parameters and initial condition|
|`observalbe.jl`|Model observables for correlating simulation results with experimental observations|
|`simulation.jl`|Simulation condition|
|`experimental_data.jl`|Experimental measurements for determining parameters|
|`set_search_param.jl`|Model parameters to optimize and search region|
|`fitness.jl`|An objective function to be minimized, i.e., the distance between model simulation and experimental data|

**optimize**(```model```::ExecModel, ```index_of_parameter_set```::Int; ```max_generation```::Int=10000, ```allowable_error```::Float64=0.0)
- Finds a parameter set that reproduces experimental observations.

## Estimate unknown model parameters

```julia
using BioMASS

model = load_model("./fos_model")

optimize(model, 1, max_generation=20000, allowable_error=0.5)
```
## Simultaneous parameter optimization

### Using module ```Distributed```
```julia
using Distributed
addprocs(); # add worker processes
@everywhere using BioMASS

@everywhere begin
    model = load_model("./fos_model")
    function optimize_parallel(i)
        optimize(model, i, max_generation=20000, allowable_error=0.5)
    end
end

pmap(optimize_parallel, 1:10)
```

### Calling multiple bash scripts

- main.jl

```julia
using BioMASS

model = load_model("./fos_model")

if abspath(PROGRAM_FILE) == @__FILE__
    optimize(model, 1, max_generation=20000, allowable_error=0.5)
end
```

- optimize_parallel.sh

```bash
#!/bin/sh

for i in $(seq 1 10); do
    nohup julia main.jl $i >> errout/$i.log  2>&1 &
done

# To terminate the process,
# $ pgrep -f main.jl | xargs kill -9
```

Run optimize_parallel.sh

```bash
$ mkdir errout
$ sh optimize_parallel.sh
```