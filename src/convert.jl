using PyCall

function __init__()
    py"""
    import os
    import shutil
    import re


    def main(model_path):
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
        try:
            import numpy as np
        except ImportError:
            print("numpy: Not installed")
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
    """
end

param2biomass(model_path::String) = py"main"(model_path)