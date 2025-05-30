module DataLoads

export load_data, StructAllData

# load functions
using ..DataLoadsFunc, ..DataAdjustments

# load functions from packages
import CSV: CSV
import DataFrames: DataFrame
import ..Params: StructAllParams

mutable struct StructAllData
    RWParams::StructRWParams
    regionParams::StructRWParams
    FFsupplyCurves::StructFFsupplyCurves 
    GsupplyCurves::StructGsupply
    projectionswind::Matrix{Float64}
    projectionssolar::Matrix{Float64}
    curtmat::Array{Float64,3}
    batteryrequirements::Matrix{Float64}
    sectoralempshares::Matrix{Union{Float64, Missing}}
    samplepointssolar::Array{Float64, 3}
    samplepointswind::Array{Float64, 3}
    samplepointsbat::Array{Float64, 3}
    R_LR::Float64
    wage_init::Vector{Float64}
    KR_init_S::Matrix{Float64}
    KR_init_W::Matrix{Float64}
end


function load_data(P::StructAllParams, D::String)
    # initialize data
    wage_init = Vector{Float64}(undef, 2531)
    secshares = Matrix{Float64}(undef, 10, 2)
    sectoralempshares = Matrix{Union{Float64, Missing}}(undef, 2531, 11)
    curtmatmx = Matrix{Float64}(undef, 21, 21)
    curtmat = zeros(21, 21, 3)
    batteryrequirements = Matrix{Float64}(undef, 15, 2)


    # ---------------------------------------------------------------------------- #
    #                                   Load Data                                  #
    # ---------------------------------------------------------------------------- #

    # load global csv data

    #regions_all, Linedata, majorregions_all, Fossilcapital, Renewablecapital, Electricityregionprices = load_csv_data(D)
    # variables not used later: regions_all, Fossilcapital, Renewablecapital, Electricityregionprices

    Linedata = CSV.File("$D/ModelDataset/Linedata.csv") |> DataFrame
    majorregions_all = CSV.File("$D/ModelDataset/majorregions.csv") |> DataFrame



    # initiate wage
    w_i!(wage_init, P.regions)

    # long run interest rate
    R_LR = 1/P.params.beta

    rP_LR = R_LR - 1 + P.params.deltaB    # long run return on production capital
    # this variable is calculated in SteadyState.jl differently. This one is only used locally to this file

    # create RWParams
    RWParams = fill_RWParams(majorregions_all, P.majorregions, P.regions, Linedata, P.params, wage_init, P.thetaS, P.thetaW, P.popelas, rP_LR, D);

    # load in sectoral shares

    sec_shares!(secshares, sectoralempshares, D)

    # ---------------------------------------------------------------------------- #
    #                             Load Fossil Fuel Data                            #
    # ---------------------------------------------------------------------------- #

    # FFsupplyCurves
    FFsupplyCurves = fill_FFsupply(D)

    # GsupplyCurves
    GsupplyCurves = fill_Gsupply(D)

    P.regions.reserves .= P.regions.reserves ./ sum(P.regions.reserves)

    # ---------------------------------------------------------------------------- #
    #                            Convert to Model Inputs                           #
    # ---------------------------------------------------------------------------- #

    # regionparams
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

    # linear map
    # storagey, storagex = linear_map()
    # have not seen either variable called anwhere else

    # create curtmat

    curtmat, samplepointssolar, samplepointswind, samplepointsbat = create_curtmat!(curtmatmx, curtmatmx, curtmatmx, curtmat, D);

    # import battery requirements
    battery_req!(batteryrequirements, D)

    return StructAllData(
        RWParams,
        regionParams,
        FFsupplyCurves,
        GsupplyCurves,
        projectionswind,
        projectionssolar,
        curtmat,
        batteryrequirements,
        sectoralempshares,
        samplepointssolar,
        samplepointswind,
        samplepointsbat,
        R_LR,
        wage_init,
        KR_init_S,
        KR_init_W
    )
end

end