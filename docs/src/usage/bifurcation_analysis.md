# Bifurcation Analysis

A numerical study of the changes in the dynamics and stability of a system upon variations in its parameters.

## Procedure for stability analysis at fixed points

Consider the following system of ordinary differential equations: d***x***/dt = ***F***(***x***)
1. Determine the fixed point vector, ***x***<sup>∗</sup>, solving ***F***(***x***<sup>∗</sup>) = **0**
1. Construct the Jacobian matrix, __*J*__(***x***) = *∂*__*F*__(***x***)/*∂*__*x*__
1. Compute eigenvalues of __*J*__(***x***<sup>∗</sup>): det |__*J*__(***x***<sup>∗</sup>) − *λ*__*E*__| = 0
1. Conclude on stability or instability of ***x***<sup>∗</sup> based on the real parts of eigenvalues
    - All eigenvalues have real parts less than zero → ***x***<sup>∗</sup> is stable
    - At least one of the eigenvalues has a real part greater than zero → ***x***<sup>∗</sup> is unstable

## Usage
See [examples/bifurcation](https://github.com/himoto/BioMASS.jl/tree/master/examples/bifurcation).