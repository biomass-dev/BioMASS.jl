using PyCall

py"""
import os
import shutil
import re


def main(model_path):
    '''Convert fitparam/n/*.dat -> out/n/*.npy

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
        print('numpy: Not installed')
    
    n_file = []
    fitparam_files = os.listdir(model_path.strip('/') + '/fitparam')
    for file in fitparam_files:
        if re.match(r'\d', file):
            n_file.append(int(file))
    for nth_paramset in n_file:
        os.makedirs(
            model_path.strip('/') 
            + '/dat2npy/out/{:d}'.format(nth_paramset), exist_ok=True
        )
        nth_fitparam_files = os.listdir(
            model_path.strip('/') + '/fitparam/{:d}'.format(nth_paramset)
        )
        for dat_file in nth_fitparam_files:
            if 'fit' in dat_file:
                '''
                - fit_param%d.dat -> fit_param%d.npy
                - best_fitness.dat -> best_fitness.npy
                '''
                try:
                    data = np.loadtxt(
                        model_path.strip('/') + '/fitparam/{:d}/{}'.format(
                            nth_paramset, dat_file
                        ), dtype='float'
                    )
                except ValueError:
                    pass
            else:
                '''
                - count_num.dat -> count_num.npy
                - generation.dat -> generation.npy
                '''
                data = np.loadtxt(
                    model_path.strip('/') + '/fitparam/{:d}/{}'.format(
                        nth_paramset, dat_file
                    ), dtype='int'
                )
            np.save(
                model_path.strip('/') + '/dat2npy/out/{:d}/'.format(nth_paramset)
                + dat_file.replace('.dat', '.npy'), data
            )
        if os.path.isfile(
                model_path.strip('/') + '/logs/{:d}.log'.format(nth_paramset)):
            shutil.copyfile(
                model_path.strip('/') + '/logs/{:d}.log'.format(nth_paramset),
                model_path.strip('/') 
                + '/dat2npy/out/{:d}/optimization.log'.format(nth_paramset)
            )
"""

function param2biomass(model_path::String)
    try
        py"main"(model_path)
    catch
        return
    end
end