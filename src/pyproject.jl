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


    def convert(model_path):
        ''' Convert fitparam/n/*.dat -> out/n/*.npy

        Parameters
        ----------
        model_path : str
            Path to your model written in Julia.

        Usage
        -----
        $ python dat2npy.py
        $ mv dat2npy/out/ path_to_biomass/

        '''
        n_file = []
        fitparam_files = os.listdir(
            os.path.join(
                model_path,
                "fitparam",
            )
        )
        for file in fitparam_files:
            if re.match(r"\d", file):
                n_file.append(int(file))
        for nth_paramset in n_file:
            os.makedirs(
                os.path.join(
                    model_path,
                    "dat2npy",
                    "out",
                    f"{nth_paramset:d}",
                ),
                exist_ok=True,
            )
            nth_fitparam_files = os.listdir(
                os.path.join(
                    model_path,
                    "fitparam",
                    f"{nth_paramset:d}",
                )
            )
            for dat_file in nth_fitparam_files:
                if os.path.splitext(dat_file)[-1] == ".dat":
                    if "fit" in dat_file:
                        '''
                        - fit_param%d.dat -> fit_param%d.npy
                        - best_fitness.dat -> best_fitness.npy
                        '''
                        try:
                            data = np.loadtxt(
                                os.path.join(
                                    model_path,
                                    "fitparam",
                                    f"{nth_paramset:d}",
                                    f"{dat_file}",
                                ),
                                dtype="float",
                            )
                        except ValueError:
                            pass
                    else:
                        '''
                        - count_num.dat -> count_num.npy
                        - generation.dat -> generation.npy
                        '''
                        data = np.loadtxt(
                            os.path.join(
                                model_path,
                                "fitparam",
                                f"{nth_paramset:d}",
                                f"{dat_file}",
                            ),
                            dtype="int",
                        )
                np.save(
                    os.path.join(
                        model_path,
                        "dat2npy",
                        "out",
                        f"{nth_paramset:d}",
                        dat_file.replace(".dat", ".npy"),
                    ),
                    data,
                )
            if os.path.isfile(
                os.path.join(
                    model_path,
                    "logs",
                    f"{nth_paramset:d}.log",
                )
            ):
                shutil.copyfile(
                    os.path.join(
                        model_path,
                        "logs",
                        f"{nth_paramset:d}.log",
                    ),
                    os.path.join(
                        model_path,
                        "dat2npy",
                        "out",
                        f"{nth_paramset:d}",
                        "optimization.log",
                    ),
                )


    class _Logger(object):

        def __init__(self, model_path: str, x_id: int):
            self.log = open(
                os.path.join(model_path, "out", DIRNAME + str(x_id), "optimization.log"),
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

            self.savedir = os.path.join(self.model_path, "out", f"{self.x_id}")
            if os.path.isdir(self.savedir):
                raise ValueError(
                    f"out{os.sep}{self.x_id} already exists in {self.model_path}. "
                    "Use another parameter id."
                )
            else:
                os.makedirs(self.savedir)
            os.makedirs(os.path.join(self.model_path, "out", DIRNAME + str(self.x_id)), exist_ok=True)
            self.default_stdout = sys.stdout

        def minimize(self, *args, **kwargs):
            try:
                from scipy.optimize import differential_evolution
            except ImportError:
                print("scipy: Not installed.")
            os.makedirs(os.path.join(self.model_path, "out", DIRNAME + str(self.x_id)), exist_ok=True)
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
                os.path.join(self.model_path, "out", DIRNAME + str(self.x_id), "optimization.log"),
                self.savedir,
            )

            best_fitness: float = self.model_objective(x)
            n_iter = self._get_n_iter()
            np.save(os.path.join(self.savedir, "best_fitness"), best_fitness)
            np.save(os.path.join(self.savedir, "count_num"), n_iter)
            np.save(os.path.join(self.savedir, "generation"), n_iter)
            np.save(os.path.join(self.savedir, f"fit_param{n_iter}"), x)
            if cleanup:
                shutil.rmtree(os.path.join(self.model_path, "out", DIRNAME + str(self.x_id)))


    def optimize(
        model_path,
        model_objective,
        model_gene2val,
        n_search_param,
        x_id,
        *,
        optimizer_options: Optional[dict] = None,
    ) -> None:
        if optimizer_options is None:
            optimizer_options = {}
        optimizer_options.setdefault("strategy", "best1bin")
        optimizer_options.setdefault("maxiter", 50)
        optimizer_options.setdefault("popsize", 3)
        optimizer_options.setdefault("tol", 1e-4)
        optimizer_options.setdefault("mutation", 0.1)
        optimizer_options.setdefault("recombination", 0.5)
        optimizer_options.setdefault("disp", True)
        optimizer_options.setdefault("polish", False)
        optimizer_options.setdefault("workers", 1)

        if not optimizer_options["disp"]:
            raise ValueError("Set optimizer_options['disp'] to True.")

        optimizer = Optimizer(
            model_path,
            model_objective,
            model_gene2val,
            x_id,
        )
        res = optimizer.minimize(
            model_objective,
            [(0.0, 1.0) for _ in range(n_search_param)],
            **optimizer_options,
        )
        param_values = model_gene2val(res.x)
        optimizer.import_solution(param_values)
    """
end


param2biomass(model_path::String) = py"convert"(model_path)


function scipy_differential_evolution(
    model::Model,
    x_id::Int,
    optimizer_options::Union{Dict,Nothing}=nothing
)::Nothing
    search_bounds::Matrix{Float64} = model.search_region()
    n_search_param::Int = size(search_bounds)[2]
    return py"optimize"(
        model.path, model.obj_func, model.gene2val, n_search_param, x_id,
        optimizer_options=optimizer_options
    )
end