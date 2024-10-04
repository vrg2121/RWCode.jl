module Params
export setup_parameters
using ..ParamsFunctions

"""
    setup_parameters(Data::String, Guesses::String)

Set parameters for the model run.

## Inputs
- `Data` -- full path to Data folder.
- `Guesses` -- full path to Guesses folder.

## Output
Model parameters wrapped in a named tuple. The named tuple 
    contains constants and structs like `params`.

## Notes
Descriptions for each parameter are in the submodule `ParamsFunctions.jl` comments.

## Example
```julia-repl
julia> D = "path/to/Data"
julia> G = "path/to/Guesses"
julia> P = setup_parameters(D, G) # the named tuple is saved in variable P
```
"""
function setup_parameters(D::String, G::String)

    thetaS = Vector{Float64}(undef, 2531)
    thetaW = Vector{Float64}(undef, 2531)

    regions, majorregions, Linecounts = load_parameters_csv(D)

    params = fill_params(regions, majorregions, Linecounts, D, G)

    capacityfactorS = 1

    thetaS = tS!(thetaS, capacityfactorS, regions)

    capacityfactorW = 1

    thetaW = tW!(thetaW, capacityfactorW, regions)
    theta = max.(thetaS, thetaW)

    # update weights for iterations
    updw_w = 0.1
    upw_z = 0.001
    updwk = 0.5
    updwF = 0.001
    upda = 0.001
    #updw_D = 0.01
    #updw_l = 0.00001
    #updw_k = 0.01
    #upw_zmax = 0.001
    #updw_alpha = 0.2
    #updps = 0.1

    # adjust for population costs
    popelas = 0.2

    # speed up for line constraints to ignore
    linconscount = 200

    # relative advantage of fossil fuels initially
    kappa = 1

    # exogenous growth rate in extraction productivity
    g = 0.01

    # variable initiliazation 
    pkwh_B = 412.37         # Initial price of batteries
    pkw_solar = 2000        # Initial price of solar
    pB_shifter = 0.0    
    curtailmentswitch = 0   # Initial curtailment

    # Transition Parameters
    T = 500
    decayp = -0.05

    return (
        params = params,
        thetaS = thetaS,
        theta = theta,
        thetaW = thetaW,
        regions = regions,
        majorregions = majorregions,
        popelas = popelas,
        T = T,
        Linecounts = Linecounts,
        linconscount = linconscount,
        kappa = kappa,
        updw_w = updw_w,
        upw_z = upw_z,
        curtailmentswitch = curtailmentswitch,
        decayp = decayp,
        pB_shifter = pB_shifter,
        pkw_solar = pkw_solar,
        pkwh_B = pkwh_B,
        g = g,
        upda = upda,
        updwF = updwF,
        updwk = updwk
    )

end
end