## Parameter Estimation

#### **load_model**(```path_to_model```::String)
- Loads a BioMASS model. The model must include the following files:

|Name|Content|
|---|---|
|`name2idx/`|Names of model parameters and species|
|`set_model.jl`|Differential equation, parameters and initial condition|
|`observalbe.jl`|Model observables|
|`simulation.jl`|Simulation condition|
|`experimental_data.jl`|Experimental data|
|`set_search_param.jl`|Model parameters to optimize and search region|
|`fitness.jl`|An objective function to be minimized, i.e., the distance between model simulation and experimental data|


#### **optimize**(```model```::ExecModel, ```index_of_parameter_set```::Int, ```max_generation```::Int, ```allowable_error```::Float64)
- Finds a parameter set that reproduces experimental observations.

### Estimate unknown model parameters against experimental observations
```julia
using BioMASS

model = load_model("./fos_model")

optimize(model, 1, max_generation=20000, allowable_error=0.5)
```
If you want to search multiple parameter sets (from 1 to 10) simultaneously,
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