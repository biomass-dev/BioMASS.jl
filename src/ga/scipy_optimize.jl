module SciPyOptimize

export fmin_powell

using PyCall

function __init__()
    py"""
    import warnings
    import numpy as np
    from scipy.optimize import minimize
    from typing import Callable


    warnings.filterwarnings('ignore')
    
    def modified_powell(
            objective: Callable,
            n_gene: int,
            population: np.ndarray,
            ip: np.ndarray
    ) -> np.ndarray:
        lower = np.min(population[ip, :n_gene], axis=0)
        upper = np.max(population[ip, :n_gene], axis=0)
        direc = np.identity(n_gene) * 0.3 * (upper - lower)
        res = minimize(
            objective,
            population[ip[0], :n_gene],
            method='Powell',
            bounds=tuple(zip(lower, upper)),
            options={
                'maxiter' : 10,
                'maxfev' : 100,
                'direc' : direc,
            }
        )
        population[ip[0], :n_gene] = res.x
        population[ip[0], -1] = objective(res.x)

        return population
    """
end

function fmin_powell(
        objective::Function,
        n_gene::Int,
        population::Matrix{Float64},
        ip::Vector{Int}
)::Matrix{Float64}
    return py"modified_powell"(objective, n_gene, population, ip.-1)
end

end # module