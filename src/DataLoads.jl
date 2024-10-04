module DataLoads

export load_data

using ..DataLoadsFunc, ..DataAdjustments

import CSV: CSV
import DataFrames: DataFrame

"""
   load_data(P::NamedTuple, D::String)

Loads and configures all data for the model
   
## Inputs
- `P::NamedTuple` -- NamedTuple of model parameters. Output from `P = setup_parameters(Data, Guesses)`
- `D::String` -- path to Data folder. `D = "path/to/Data"`

## Outputs
NamedTuple of all data necessary to run the model. The NamedTuple contains constants 
    and structs like RWParams and regionParams.

## Notes
Descriptions of the data are in the submodule `DataLoadsFunc.jl` comments.

## Example
```julia-repl
julia> D = "path/to/Data"
julia> G = "path/to/Guesses"
julia> P = setup_parameters(D, G)
julia> DL = load_data(P, D) # all data is saved in named tuple variable, DL

```
"""
function load_data(P::NamedTuple, D::String)
    wage_init = Vector{Float64}(undef, 2531)
    secshares = Matrix{Float64}(undef, 10, 2)
    sectoralempshares = Matrix{Union{Float64, Missing}}(undef, 2531, 11)
    curtmatmx = Matrix{Float64}(undef, 21, 21)
    curtmat = zeros(21, 21, 3)
    batteryrequirements = Matrix{Float64}(undef, 15, 2)

    regions_all, Linedata, majorregions_all, Fossilcapital, Renewablecapital, Electricityregionprices = load_csv_data(D)

    wage_init = w_i!(wage_init, P.regions)

    # long run interest rate
    R_LR = 1/P.params.beta

    rP_LR = R_LR - 1 + P.params.deltaB    # long run return on production capital
    # this variable is calculated in SteadyState.jl differently. This one is only used locally to this file

    RWParams = fill_RWParams(majorregions_all, P.majorregions, P.regions, Linedata, P.params, wage_init, P.thetaS, P.thetaW, P.popelas, rP_LR, D);

    secshares, sectoralempshares = sec_shares!(secshares, sectoralempshares, D)

    # obtain country supply curves
    FFsupplyCurves = fill_FFsupply(D)
    GsupplyCurves = fill_Gsupply(D)

    P.regions.reserves .= P.regions.reserves ./ sum(P.regions.reserves)

    regionParams = deepcopy(RWParams)

    # initial renewable capital
    KR_init_S = RWParams.KR ./ 2
    KR_init_W = RWParams.KR ./ 2

    # ---------------------------------------------------------------------------- #
    #                        LOAD EXOGENOUS TECH PROJECTIONS                       #
    # ---------------------------------------------------------------------------- #

    projectionswind = DataFrame(CSV.File("$D/ModelDataset/Windprojectionsdata.csv")) |> Matrix
    projectionssolar = DataFrame(CSV.File("$D/ModelDataset/Solarprojectionsdata.csv")) |> Matrix

    # ---------------------------------------------------------------------------- #
    #                             LOAD CURTAILMENT DATA                            #
    # ---------------------------------------------------------------------------- #

    # linear map (unused variables)
    #storagey, storagex = linear_map()

    curtmat, samplepointssolar, samplepointswind, samplepointsbat = create_curtmat!(curtmatmx, curtmatmx, curtmatmx, curtmat, D);

    # import battery requirements
    batteryrequirements = battery_req!(batteryrequirements, D)

    return(
        RWParams = RWParams,
        regionParams = regionParams,
        FFsupplyCurves = FFsupplyCurves,
        GsupplyCurves = GsupplyCurves,
        projectionswind = projectionswind,
        projectionssolar = projectionssolar,
        curtmat = curtmat,
        batteryrequirements = batteryrequirements,
        sectoralempshares = sectoralempshares,
        samplepointssolar = samplepointssolar,
        samplepointswind = samplepointswind,
        samplepointsbat = samplepointsbat,
        R_LR = R_LR, 
        wage_init = wage_init, 
        KR_init_S = KR_init_S, 
        KR_init_W = KR_init_W
        )
end
end