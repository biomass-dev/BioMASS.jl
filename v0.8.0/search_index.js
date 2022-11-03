var documenterSearchIndex = {"docs":
[{"location":"references/","page":"References","title":"References","text":"Imoto, H., Zhang, S. & Okada, M. A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data—Application to the ErbB Receptor Signaling Pathway. Cancers 12, 2878 (2020). https://doi.org/10.3390/cancers12102878\nImoto, H., Yamashiro, S. & Okada, M. A text-based computational framework for patient -specific modeling for classification of cancers. iScience 25, 103944 (2022). https://doi.org/10.1016/j.isci.2022.103944\nNakakuki, T. et al. Ligand-specific c-Fos expression emerges from the spatiotemporal control of ErbB network dynamics. Cell 141, 884–896 (2010). https://doi.org/10.1016/j.cell.2010.03.054\nYao, G., Lee, T. J., Mori, S., Nevins, J. R. & You, L. A bistable Rb-E2F switch underlies the restriction point. Nat. Cell Biol. 10, 476–482 (2008). https://doi.org/10.1038/ncb1711\nBarr, A. R., Heldt, F. S., Zhang, T., Bakal, C. & Novák, B. A Dynamical Framework for the All-or-None G1/S Transition. Cell Syst. 2, 27–37 (2016). https://doi.org/10.1016/j.cels.2016.01.001\nRata, S. et al. Two Interlinked Bistable Switches Govern Mitotic Control in Mammalian Cells. Curr. Biol. 28, 3824-3832.e6 (2018). https://doi.org/10.1016/j.cub.2018.09.059","category":"page"},{"location":"usage/model_construction/#Model-Construction","page":"Model Construction","title":"Model Construction","text":"","category":"section"},{"location":"usage/model_construction/","page":"Model Construction","title":"Model Construction","text":"pasmopy.Text2Model allows you to build a BioMASS model from text. You simply describe biochemical reactions and the molecular mechanisms extracted from text are converted into an executable model.","category":"page"},{"location":"usage/model_construction/#Example","page":"Model Construction","title":"Example","text":"","category":"section"},{"location":"usage/model_construction/","page":"Model Construction","title":"Model Construction","text":"This example shows you how to build a simple Michaelis-Menten two-step enzyme catalysis model with Pasmopy.","category":"page"},{"location":"usage/model_construction/","page":"Model Construction","title":"Model Construction","text":"E + S ⇄ ES → E + P","category":"page"},{"location":"usage/model_construction/","page":"Model Construction","title":"Model Construction","text":"An enzyme, E, binding to a substrate, S, to form a complex, ES, which in turn releases a product, P, regenerating the original enzyme.","category":"page"},{"location":"usage/model_construction/","page":"Model Construction","title":"Model Construction","text":"Prepare a text file describing biochemical reactions (e.g., michaelis_menten.txt)\n```  E binds S <–> ES | kf=0.003, kr=0.001 | E=100, S=50  ES dissociates to E and P | kf=0.002, kr=0\n@obs Substrate: u[S]  @obs Efree: u[E]  @obs Etotal: u[E] + u[ES]  @obs Product: u[P]  @obs Complex: u[ES]\n@sim tspan: [0, 100]  ```\nConvert the text into an executable model\nshell  $ python  ```python\nfrom pasmopy import Text2Model description = Text2Model(\"michaelismenten.txt\", lang=\"julia\") description.convert()  # generate 'michaelismenten_jl/'\nModel information\n2 reactions  4 species  4 parameters  ```\nRun simulation\nshell  $ julia  ```julia  using BioMASS\nmodel = Model(\"./michaelismentenjl\");  run_simulation(model)  ```\n(Image: michaelis_menten)","category":"page"},{"location":"usage/bifurcation_analysis/#Bifurcation-Analysis","page":"Bifurcation Analysis","title":"Bifurcation Analysis","text":"","category":"section"},{"location":"usage/bifurcation_analysis/","page":"Bifurcation Analysis","title":"Bifurcation Analysis","text":"A numerical study of the changes in the dynamics and stability of a system upon variations in its parameters.","category":"page"},{"location":"usage/bifurcation_analysis/","page":"Bifurcation Analysis","title":"Bifurcation Analysis","text":"(Image: )","category":"page"},{"location":"usage/bifurcation_analysis/#Procedure-for-stability-analysis-at-fixed-points","page":"Bifurcation Analysis","title":"Procedure for stability analysis at fixed points","text":"","category":"section"},{"location":"usage/bifurcation_analysis/","page":"Bifurcation Analysis","title":"Bifurcation Analysis","text":"Consider the following system of ordinary differential equations:","category":"page"},{"location":"usage/bifurcation_analysis/","page":"Bifurcation Analysis","title":"Bifurcation Analysis","text":"dfracdxdt = F(x)","category":"page"},{"location":"usage/bifurcation_analysis/","page":"Bifurcation Analysis","title":"Bifurcation Analysis","text":"Determine the fixed point vector, x^*, solving F(x^*) = 0\nConstruct the Jacobian matrix, J(x) = dfracpartial F(x)partial x\nCompute eigenvalues of J(x^*): J(x^*)  λE = 0\nConclude on stability or instability of x^* based on the real parts of eigenvalues\nAll eigenvalues have real parts less than zero → x^* is stable\nAt least one of the eigenvalues has a real part greater than zero → x^* is unstable","category":"page"},{"location":"usage/bifurcation_analysis/#Usage","page":"Bifurcation Analysis","title":"Usage","text":"","category":"section"},{"location":"usage/bifurcation_analysis/","page":"Bifurcation Analysis","title":"Bifurcation Analysis","text":"See examples/bifurcation.","category":"page"},{"location":"usage/parameter_estimation/#Parameter-Estimation","page":"Parameter Estimation","title":"Parameter Estimation","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"(Image: )","category":"page"},{"location":"usage/parameter_estimation/#Core-functions","page":"Parameter Estimation","title":"Core functions","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"","category":"page"},{"location":"usage/parameter_estimation/#Model(path_to_model::String)","page":"Parameter Estimation","title":"Model(path_to_model::String)","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Load a BioMASS model. The model must include the following files:","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Name Content\nname2idx/ Names of model parameters and species\node.jl Differential equation, parameters and initial condition\nobservalbe.jl Model observables for correlating simulation results with experimental observations\nsimulation.jl Simulation condition\nexperimental_data.jl Experimental measurements\nsearch_param.jl Lower and upper bounds of model parameters to be estimated\nproblem.jl An objective function to be minimized, i.e., the distance between model simulation and experimental data","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Parameters\npath_to_model::String\nThe model folder to read.\nReturns\nmodel::Model\nThe executable model in BioMASS.","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"note: Note\npasmopy.Text2Model allows you to build a BioMASS model from text [Imoto et al., 2022]. You simply describe biochemical reactions and the molecular mechanisms extracted from text are converted into an executable model. To build a model for BioMASS.jl, please set lang=\"julia\".","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"","category":"page"},{"location":"usage/parameter_estimation/#scipy_differential_evolution(model::Model,-ix_id::Int,-kwargs...)","page":"Parameter Estimation","title":"scipy_differential_evolution(model::Model, ix_id::Int, kwargs...)","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Estimate model parameters from experimental data.","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Parameters\nmodel::Model\nModel for parameter estimation.\nx_id::Int\nIndex of parameter set to estimate.\nkwargs...\nKeyword arguments to pass to scipy.optimize.differential_evolution.","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"","category":"page"},{"location":"usage/parameter_estimation/#run_simulation(model::Model,-viz_type::String,-show_all::Boolfalse,-stdev::Boolfalse)","page":"Parameter Estimation","title":"run_simulation(model::Model, viz_type::String, show_all::Bool=false, stdev::Bool=false)","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Save simulation results with optimized parameter values.","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Parameters\nviz_type::String\n\"average\"\n\"best\"\n\"original\"\n\"experiment\"\nshow_all::Bool (default: false)\nWhether to show all simulation results.\nstdev::Bool (default: false)\nIf True, the standard deviation of simulated values will be shown (only available for \"average\" visualization type).\nsave_format::String (default: \"pdf\")\nEither \"png\" or \"pdf\", indicating whether to save figures as png or pdf format.","category":"page"},{"location":"usage/parameter_estimation/#Estimate-unknown-model-parameters","page":"Parameter Estimation","title":"Estimate unknown model parameters","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"using BioMASS\n\nmodel = Model(\"./examples/fos_model\");\ninitpop = generate_initial_population(model)\nscipy_differential_evolution(model, 1, init=initpop)","category":"page"},{"location":"usage/parameter_estimation/#Simultaneous-parameter-optimization","page":"Parameter Estimation","title":"Simultaneous parameter optimization","text":"","category":"section"},{"location":"usage/parameter_estimation/#Using-module-Distributed","page":"Parameter Estimation","title":"Using module Distributed","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"using Distributed\naddprocs(); # add worker processes\n@everywhere using BioMASS\n\n@everywhere begin\n    model = Model(\"./examples/fos_model\")\n    function optimize_parallel(i)\n        scipy_differential_evolution(model, i)\n    end\nend\n\npmap(optimize_parallel, 1:10)","category":"page"},{"location":"usage/parameter_estimation/#Calling-multiple-bash-scripts","page":"Parameter Estimation","title":"Calling multiple bash scripts","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"main.jl","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"using BioMASS\n\nmodel = Model(\"./examples/fos_model\")\n\nif abspath(PROGRAM_FILE) == @__FILE__\n    scipy_differential_evolution(model, parse(Int64, ARGS[1]))\nend","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"optimize_parallel.sh","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"#!/bin/sh\n\nfor i in $(seq 1 10); do\n    nohup julia main.jl $i >> errout/$i.log  2>&1 &\ndone\n\n# To terminate the process,\n# $ pgrep -f main.jl | xargs kill -9","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"Run optimize_parallel.sh","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"$ mkdir errout\n$ sh optimize_parallel.sh","category":"page"},{"location":"usage/parameter_estimation/#How-to-track-optimization-process","page":"Parameter Estimation","title":"How to track optimization process","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"The temporary result will be saved in path_to_model/fitparam/n/optimization.log.","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"$ tail examples/fos_model/fitparam/1/optimization.log","category":"page"},{"location":"usage/parameter_estimation/#Visualization-of-simulation-results","page":"Parameter Estimation","title":"Visualization of simulation results","text":"","category":"section"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"The simulation results will be saved in figure/.","category":"page"},{"location":"usage/parameter_estimation/","page":"Parameter Estimation","title":"Parameter Estimation","text":"run_simulation(model, viz_type=\"best\", show_all=true)","category":"page"},{"location":"#BioMASS.jl","page":"Home","title":"BioMASS.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: License: MIT) (Image: Source) (Image: Actions Status) (Image: Cancers Paper)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This module provides a Julia interface to the BioMASS parameter estimation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The open access publication describing BioMASS is available here:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Imoto, H., Zhang, S. & Okada, M. A Computational Framework for Prediction and Analysis of Cancer Signaling Dynamics from RNA Sequencing Data—Application to the ErbB Receptor Signaling Pathway. Cancers 12, 2878 (2020). https://doi.org/10.3390/cancers12102878","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"] add BioMASS","category":"page"},{"location":"","page":"Home","title":"Home","text":"    Pages = [\n        \"usage/parameter_estimation.md\",\n        \"usage/bifurcation_analysis.md\",\n    ]\n    Depth = 3","category":"page"}]
}