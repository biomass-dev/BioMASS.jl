# Model Construction

[`pasmopy.Text2Model`](https://pasmopy.readthedocs.io/en/latest/model_development.html) allows you to build a BioMASS model from text. You simply describe biochemical reactions and the molecular mechanisms extracted from text are converted into an executable model.

## Example

This example shows you how to build a simple Michaelis-Menten two-step enzyme catalysis model with Pasmopy.

> E + S ⇄ ES → E + P

_An enzyme, E, binding to a substrate, S, to form a complex, ES, which in turn releases a product, P, regenerating the original enzyme._

1. Prepare a text file describing biochemical reactions (e.g., `michaelis_menten.txt`)

    ```
    E binds S <--> ES | kf=0.003, kr=0.001 | E=100, S=50
    ES dissociates to E and P | kf=0.002, kr=0

    @obs Substrate: u[S]
    @obs E_free: u[E]
    @obs E_total: u[E] + u[ES]
    @obs Product: u[P]
    @obs Complex: u[ES]

    @sim tspan: [0, 100]
    ```

1. Convert the text into an executable model

    ```shell
    $ python
    ```
    ```python
    >>> from pasmopy import Text2Model
    >>> description = Text2Model("michaelis_menten.txt", lang="julia")
    >>> description.convert()  # generate 'michaelis_menten_jl/'
    Model information
    -----------------
    2 reactions
    4 species
    4 parameters
    ```

1. Run simulation

    ```shell
    $ julia
    ```
    ```julia
    using BioMASS

    model = Model("./michaelis_menten_jl");
    run_simulation(model)
    ```

    ![michaelis_menten](https://raw.githubusercontent.com/pasmopy/pasmopy/master/docs/_static/img/michaelis_menten_sim.png)