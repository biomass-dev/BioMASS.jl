using PyCall

function __init__()
    py"""
    import os
    import shutil
    import re
    import sys
    import warnings
    from typing import Callable, Optional, List, Union

    try:
        import numpy as np
    except ImportError:
        print("numpy: Not installed")


    DIRNAME = "_tmp"


    class _Logger(object):

        def __init__(self, model_path: str, x_id: int):
            self.log = open(
                os.path.join(model_path, "fitparam", DIRNAME + str(x_id), "optimization.log"),
                mode="w",
                encoding="utf-8",
            )

        def write(self, message: str):
            self.log.write(message)


    class Optimizer(object):

        def __init__(
            self,
            model_path,
            model_objective,
            model_gene2val,
            x_id,
        ):
            self.model_path = model_path
            self.model_objective = model_objective
            self.model_gene2val = model_gene2val
            self.x_id = x_id

            self.savedir = os.path.join(self.model_path, "fitparam", f"{self.x_id}")
            if os.path.isdir(self.savedir):
                raise ValueError(
                    f"out{os.sep}{self.x_id} already exists in {self.model_path}. "
                    "Use another parameter id."
                )
            else:
                os.makedirs(self.savedir)
            os.makedirs(os.path.join(self.model_path, "fitparam", DIRNAME + str(self.x_id)), exist_ok=True)
            self.default_stdout = sys.stdout

        def minimize(self, *args, **kwargs):
            try:
                from scipy.optimize import differential_evolution
            except ImportError:
                print("scipy: Not installed.")
            os.makedirs(os.path.join(self.model_path, "fitparam", DIRNAME + str(self.x_id)), exist_ok=True)
            try:
                sys.stdout = _Logger(self.model_path, self.x_id)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    res = differential_evolution(*args, **kwargs)
                return res
            finally:
                sys.stdout = self.default_stdout

        def _get_n_iter(self) -> int:
            n_iter: int = 0
            path_to_log = os.path.join(self.savedir, "optimization.log")
            with open(path_to_log, mode="r", encoding="utf-8") as f:
                log_file = f.readlines()
            for message in log_file:
                if len(message.strip()) > 0:
                    n_iter += 1
            return n_iter

        def import_solution(self, x: Union[np.ndarray, List[float]], cleanup: bool = True) -> None:

            shutil.move(
                os.path.join(self.model_path, "fitparam", DIRNAME + str(self.x_id), "optimization.log"),
                self.savedir,
            )

            best_fitness: float = self.model_objective(x)
            n_iter = self._get_n_iter()
            np.save(os.path.join(self.savedir, "best_fitness"), best_fitness)
            np.save(os.path.join(self.savedir, "count_num"), n_iter)
            np.save(os.path.join(self.savedir, "generation"), n_iter)
            np.save(os.path.join(self.savedir, f"fit_param{n_iter}"), x)
            if cleanup:
                shutil.rmtree(os.path.join(self.model_path, "fitparam", DIRNAME + str(self.x_id)))


    def optimize(
        model_path,
        model_objective,
        model_gene2val,
        n_search_param,
        x_id,
        **kwargs,
    ) -> None:
        optimizer = Optimizer(
            model_path,
            model_objective,
            model_gene2val,
            x_id,
        )
        res = optimizer.minimize(
            model_objective,
            [(0.0, 1.0) for _ in range(n_search_param)],
            **kwargs,
        )
        param_values = model_gene2val(res.x)
        optimizer.import_solution(param_values)
    """
end


param2biomass(model_path::String) = py"convert"(model_path)


numpy_load(path::String) = py"np.load"(path)


function scipy_differential_evolution(
    model::Model,
    x_id::Int;
    strategy::String="best1bin",
    maxiter::Int=100,
    popsize::Int=3,
    tol::Float64=1e-4,
    mutation::Union{Float64,Tuple{Float64,Float64}}=0.1,
    recombination::Float64=0.5,
    seed::Union{Nothing,Int}=nothing,
    disp::Bool=true,
    polish::Bool=false,
    init::Union{String,Matrix{Float64}}="latinhypercube",
    atol::Float64=0.0,
    updating::String="immediate"
)::Nothing
    search_bounds::Matrix{Float64} = model.search_region()
    n_search_param::Int = size(search_bounds)[2]
    if !disp
        error("Set 'disp' to true.")
    end
    if polish
        error("Set 'polish' to false.")
    end
    return py"optimize"(
        model.path,
        model.obj_func,
        model.gene2val,
        n_search_param,
        x_id,
        strategy=strategy,
        maxiter=maxiter,
        popsize=popsize,
        tol=tol,
        mutation=mutation,
        recombination=recombination,
        seed=seed,
        disp=disp,
        polish=polish,
        init=init,
        atol=atol,
        updating=updating,
    )
end


function generate_initial_population(
    model::Model;
    popsize::Int=3,
    threshold::Float64=1e12,
    show_progress::Bool=true
)::Matrix{Float64}
    search_bounds::Matrix{Float64} = model.search_region()
    n_gene::Int = size(search_bounds)[2]
    n_population::Int = popsize * n_gene
    population::Matrix{Float64} = fill(
        Inf, (n_population, n_gene + 1)
    )
    for i = 1:n_population
        while threshold <= population[i, end]
            for j = 1:n_gene
                population[i, j] = rand()
            end
            population[i, end] = model.obj_func(population[i, 1:n_gene])
        end
        if show_progress
            print("\r$i / $n_population")
        end
    end
    population = sortslices(population, dims=1, by=x -> x[end])
    return population[:, 1:n_gene]
end
