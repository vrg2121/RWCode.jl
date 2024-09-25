
module DataLoads

# export data
export RWParams, regionParams, FFsupplyCurves, GsupplyCurves, projectionswind, projectionssolar, curtmat, 
       batteryrequirements, sectoralempshares, samplepointssolar, samplepointswind, samplepointsbat

# export constants
export R_LR

# export variables
export wage_init, KR_init_S, KR_init_W

using CSV, DataFrames, LinearAlgebra, MAT, Main.DataAdjustments, Tables, Main.Functions
using Main.Params: regions, majorregions, params, thetaS, thetaW, popelas

# benchmarking
using Profile

# load data 
function load_csv_data()
    regions_all = CSV.File("Data/ModelDataset/Regions_sorted.csv", header=true) |> DataFrame
    Linedata = CSV.File("Data/ModelDataset/Linedata.csv") |> DataFrame
    majorregions_all = CSV.File("Data/ModelDataset/majorregions.csv") |> DataFrame
    Fossilcapital = CSV.File("Data/ModelDataset/aggcap_ffl.csv")  |> DataFrame
    Renewablecapital = CSV.File("Data/ModelDataset/aggcap_rnw.csv") |> DataFrame
    Electricityregionprices = CSV.File("Data/ModelDataset/Region_prices.csv") |> DataFrame
    return regions_all, Linedata, majorregions_all, Fossilcapital, Renewablecapital, Electricityregionprices
end

# initiate variables
function w_i(regions::DataFrame)
    wage_init = Vector{Float64}(undef, 2531)
    wage_init .= regions.wages ./ (regions.wages[1])
    coalesce(wage_init, 1.0)
    wage_init[wage_init .> 2] .= 2.0
    return wage_init
end
wage_init = w_i(regions)

# long run interest rate
R_LR = 1/params.beta
rP_LR = R_LR - 1 + params.deltaB # long run return on production capital


# create a struct to store RWParams
mutable struct StructRWParams
    Kmatrix::Matrix{Float64}
    A::Vector{Matrix{Float64}}
    Adash::Vector{Matrix{Float64}}
    O::Vector{Matrix{Float64}}
    Gam::Vector{Matrix{Float64}}
    Gam2::Vector{Matrix{Float64}}
    Gam3::Vector{Matrix{Float64}}
    B::Vector{Matrix{Float64}}
    KF::Matrix{Float64}
    maxF::Matrix{Float64}
    KR::Matrix{Float64}
    Zmax::Matrix{Float64}
    thetaS::Matrix{Float64}
    thetaW::Matrix{Float64}
    scalar::Int
    Countryshifter::Matrix{Float64}
    costshifter::Matrix{Float64}
    KP::Matrix{Float64}
    z::Matrix{Float64}
    zS::Matrix{Float64}

    # Inner Constructor
    function StructRWParams()
        s1 = [(1466, 726), (1432, 754), (59, 29), (114, 52), 
        (14, 12), (23, 14), (76, 45), (10, 10), 
        (42, 25), (157, 77), (247, 124), (738, 319)]

        s2 = [(1466, 1466), (1432, 1432), (59, 59), (114, 114), 
        (14, 14), (23, 23), (76, 76), (10, 10), 
        (42, 42), (157, 157), (247, 247), (738, 738)]

        s3 = [(726, 726), (754, 754), (29, 29), (52, 52), 
        (12, 12), (14, 14), (45, 45), (10, 10), 
        (25, 25), (77, 77), (124, 124), (319, 319)]
        
        new(
        Matrix{Float64}(undef, 4382, 2531),             # KMatrix
        [Matrix{Float64}(undef, size...) for size in s1], # A
        [Matrix{Float64}(undef, size...) for size in s1], # Adash
        [Matrix{Float64}(undef, size...) for size in s2], # O
        [Matrix{Float64}(undef, size...) for size in s1], # Gam
        [Matrix{Float64}(undef, size...) for size in s1], # Gam2
        [Matrix{Float64}(undef, size...) for size in s1], # Gam3
        [Matrix{Float64}(undef, size...) for size in s3], # B
        Matrix{Float64}(undef, 2531, 1),                  # KF
        Matrix{Float64}(undef, 2531, 1),                  # maxF
        Matrix{Float64}(undef, 2531, 1),                  # KR
        Matrix{Float64}(undef, 4382, 1),                  # Zmax
        Matrix{Float64}(undef, 2531, 1),                  # thetaS
        Matrix{Float64}(undef, 2531, 1),                  # thetaW
        0,                                                # scalar
        Matrix{Float64}(undef, 2531, 1),                  # Countryshifter
        Matrix{Float64}(undef, 2531, 1),                  # costshifter
        Matrix{Float64}(undef, 2531, 1),                  # KP
        Matrix{Float64}(undef, 2531, 1),                  # z 
        Matrix{Float64}(undef, 2531, 10)                  # zS
        )
    end
end

function create_RWParams!(RWParams::StructRWParams, majorregions::DataFrame, regions::DataFrame, params, wage_init::Vector{Float64},
                          thetaS::Vector{Float64}, thetaW::Vector{Float64}, popelas, rP_LR)

    Linedata, majorregions_all = load_csv_data()

    # create kmatrix
    Kmatrix = Matrix(CSV.File("Data/ModelDataset/Kmatrix.csv", drop=[1]) |> DataFrame)
    RWParams.Kmatrix .= Matshifter(Kmatrix)

    # create A
    Kmx = Matrix(CSV.File("Data/ModelDataset/Kmatrix_1.csv") |> DataFrame)[:, 2:majorregions_all.rowid[1]]
    RWParams.A[1] .= Matshifter(Kmx)
    RWParams.Adash[1] .= Matshifter(Kmx)

    # create Omatrix
    Omatrix = CSV.File("Data/ModelDataset/Omatrix_1.csv", drop=[1]) |> DataFrame
    Omatrix = Vector(Omatrix[!,1])
    RWParams.O[1] .= Diagonal(Omatrix)

    # create Gam
    RWParams.Gam[1] .= RWParams.O[1]^(-1) * RWParams.A[1] * inv(RWParams.A[1]' * RWParams.O[1]^(-1) * RWParams.A[1])
    RWParams.Gam2[1] .= RWParams.Gam[1]

    row_maxima = vec(maximum(RWParams.Gam[1], dims=2))
    num_repeats = majorregions_all.rowid[1] - 1
    indmax = repeat(row_maxima, 1, num_repeats)

    RWParams.Gam2[1][RWParams.Gam2[1] .< indmax] .= 0

    RWParams.Gam3[1] .= RWParams.Gam[1]
    row_min = vec(minimum(RWParams.Gam[1], dims=2))
    indmin = repeat(row_min, 1, num_repeats)
    RWParams.Gam3[1][RWParams.Gam3[1] .> indmin] .= 0

    # fill the rest of A, Adash, O in a for loop
    for jj in 2:(size(majorregions_all, 1)-1)
        #stringer = "Data/ModelDataset/Kmatrix_$(jj).csv"
        #stringer2 = "Data/ModelDataset/Omatrix_$(jj).csv"
        #Kmatrix = Matrix(CSV.File(stringer) |> DataFrame)
        Kmatrix = CSV.File("Data/ModelDataset/Kmatrix_$(jj).csv") |> Tables.matrix
        Kmatrix = Kmatrix[:, majorregions_all.rowid[jj-1] + 2 : majorregions_all.rowid[jj]]
        
        #Omatrix = Matrix(CSV.File(stringer2, drop=[1]) |> DataFrame)
        Omatrix = CSV.File("Data/ModelDataset/Omatrix_$(jj).csv", drop=[1]) |> Tables.matrix
        RWParams.A[jj] .= Matshifter(Kmatrix)
        RWParams.Adash[jj] .= Matshifter(Kmatrix)
        Omatrix = Vector(Omatrix[:,1])
        RWParams.O[jj] .= Diagonal(Omatrix)
        RWParams.Gam[jj] .= RWParams.O[jj]^(-1) * RWParams.A[jj] * inv(RWParams.A[jj]' * RWParams.O[jj]^(-1) * RWParams.A[jj])

       
        indmax =  repeat(maximum(RWParams.Gam[jj], dims=2), 1, majorregions_all.n[jj] - 1)
        RWParams.Gam2[jj] .= RWParams.Gam[jj]
        RWParams.Gam2[jj][RWParams.Gam[jj] .< indmax] .= 0
    
        
        indmin = repeat(minimum(RWParams.Gam[jj], dims=2), 1, majorregions_all.n[jj] - 1)
        RWParams.Gam3[jj] .= RWParams.Gam[jj]
        RWParams.Gam3[jj][RWParams.Gam3[jj] .> indmin] .= 0   
    end

    R = size(majorregions, 1) - 1   # regions
    I = 10                          # industries

    # create B 
    for jj in 1:size(majorregions, 1) - 1
        RWParams.B[jj] .= inv(RWParams.A[jj]' * RWParams.O[jj] * RWParams.A[jj])
    end

    # fossil fuel capital
    RWParams.KF .= regions.ffl_capacity_mw ./ regions.ffl_capacity_mw[1]
    RWParams.KF .+= 10.0^(-1)

    # set max fossil fuel use to capacity
    RWParams.maxF .= RWParams.KF

    # renewable capital
    RWParams.KR .= regions.rnw_capacity_mw ./ regions.ffl_capacity_mw[1]

    # max line flows
    RWParams.Zmax .= Linedata.lcap ./ regions.ffl_capacity_mw[1]

    # renewable potential
    RWParams.thetaS .= coalesce.(thetaS, 0.2) # set places without sunlight data to very low
    replace!(RWParams.thetaS, 0 => 0.2)
    RWParams.thetaW .= coalesce.(thetaW, 0.2) # set places without wind data to very low

    #shift potential so that urban regions don't install lots of capital, then
    #scale by relative country costshifters
    RWParams.scalar = params.renewablecostscaler
    RWParams.Countryshifter .= regions[!, :costrel]

    # define costshifter
    for kk in 1:params.N
        range = majorregions[!, :rowid2][kk]:majorregions[!, :rowid][kk]
        
        RWParams.costshifter[range] = (
            (regions[!, :pop_dens][range] .^ popelas) ./ 
            (regions[!, :pop_dens][majorregions[!, :rowid2][kk]] .^ popelas)
        ) .* RWParams.scalar .* RWParams.Countryshifter[range]
    end

    # production capital
    RWParams.KP .= regions[!, :capitalperworker]
    R_normaliser = wage_init ./ (RWParams.KP .* rP_LR)
    RWParams.KP .= RWParams.KP .* R_normaliser

    RWParams.z .= params.Z
    RWParams.zS .= repeat(params.Z, 1, 10)

end

function fill_RWParams(majorregions_all::DataFrame, majorregions::DataFrame, regions::DataFrame, Linedata::DataFrame, params, wage_init::Vector{Float64}, thetaS::Vector{Float64}, 
                        thetaW::Vector{Float64}, popelas, rP_LR)
    RWParams = StructRWParams()
    create_RWParams!(RWParams, majorregions_all, majorregions, regions, Linedata, params, wage_init, thetaS, thetaW, popelas, rP_LR)
    return RWParams
end
@time RWParams = fill_RWParams(majorregions_all, majorregions, regions, Linedata, params, wage_init, thetaS, thetaW, popelas, rP_LR)
## right now RWParams is stored in the heap and accessed by a chain of addresses. 
## Alternative method to store in stack using immutable struct is in rwparams_notes.jl 


# load in sectoral shares
function sec_shares()
    secshares = CSV.File("Data/ModelDataset/secshares.csv") |> DataFrame
    secshares = Array(secshares[:, 2:3])

    sectoralempshares = CSV.read("Data/ModelDataset/locationempshares.csv", DataFrame)
    sectoralempshares = Array(sectoralempshares)
    sectoralempshares = coalesce.(sectoralempshares, 0)
    return secshares, sectoralempshares
end

@time secshares, sectoralempshares = sec_shares()

# ---------------------------------------------------------------------------- #
#                             Load Fossil Fuel Data                            #
# ---------------------------------------------------------------------------- #


struct StructFFsupplyCurves
    Q::Matrix{Float64}
    P::Matrix{Float64}

    function StructFFsupplyCurves()
        new(
            # initiate P and Q with zeroes
            zeros(Float64, 15204, 16), 
            zeros(Float64, 15204, 16)
        )
    end

end

function create_FFsupplyCurves(FFsupplyCurves::StructFFsupplyCurves)
    countriesCurves = CSV.File("Data/FossilFuels/country_curves.csv") |> DataFrame

    # obtain country supply curves
    countries = unique(countriesCurves[!, :region_name])
    totalCountries = length(countries)
    maxPoints = 15204
    for i in 1:totalCountries
        country_name = countries[i]
        newQ = countriesCurves[countriesCurves.region_name .== country_name, :Q]
        newP = countriesCurves[countriesCurves.region_name .== country_name, :P_smooth]
        
        # skip the first row to leave 0 row to top of P and Q
        FFsupplyCurves.Q[2:1+length(newQ), i] .= newQ
        FFsupplyCurves.P[2:1+length(newP), i] .= newP
        
        # repeat the last value of newQ down the column, leave the last row as 0
        if length(newQ) < maxPoints
            FFsupplyCurves.Q[(1+length(newQ)):end-1, i] .= newQ[end]
            FFsupplyCurves.P[(1+length(newP)):end-1, i] .= newP[end]
        end
    end

end

function fill_FFsupply()
    FFsupplyCurves = StructFFsupplyCurves()
    create_FFsupplyCurves(FFsupplyCurves)
    return FFsupplyCurves
end
@profile FFsupplyCurves = fill_FFsupply()

mutable struct StructGsupply
    Q::Matrix{Float64}
    P::Matrix{Float64}

    function StructGsupply()
        new(
            Matrix{Float64}(undef, 67518, 1),
            Matrix{Float64}(undef, 67518, 1)
        )
    end
end

function create_StructGsupply()

    globalCurve = CSV.File("Data/FossilFuels/global_curve.csv") |> DataFrame

    GsupplyCurves.Q[1:67517, 1] .= Vector{Float64}(globalCurve[:, :Q]) .* 100000
    GsupplyCurves.P[1:67517, 1] .= Vector{Float64}(globalCurve[:, :P_smooth]) ./ 200

    GsupplyCurves.P[67518, 1] = GsupplyCurves.P[67517, 1] * 1000
    GsupplyCurves.Q[67518, 1] = GsupplyCurves.Q[67517, 1] * 1.001

end

GsupplyCurves = StructGsupply()

create_StructGsupply()

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

projectionswind = DataFrame(CSV.File("Data/ModelDataset/Windprojectionsdata.csv")) |> Matrix
projectionssolar = DataFrame(CSV.File("Data/ModelDataset/Solarprojectionsdata.csv")) |> Matrix

# ---------------------------------------------------------------------------- #
#                             LOAD CURTAILMENT DATA                            #
# ---------------------------------------------------------------------------- #

# linear map
maxstorage=12
storepoints=100
storagey = range(0, stop=maxstorage, length=storepoints)
storagex = range(0, stop=1, length=storepoints)

# read Matrix
curtmatno = DataFrame(CSV.File("Data/CurtailmentUS/heatmap_us_mat_nostorage.csv", header = false)) |> Matrix
curtmat4 = DataFrame(CSV.File("Data/CurtailmentUS/heatmap_us_mat_4hour.csv", header = false)) |> Matrix
curtmat12 = DataFrame(CSV.File("Data/CurtailmentUS/heatmap_us_mat_12hour.csv", header = false)) |> Matrix

n = size(curtmatno, 1)
x = range(0, stop=1, length=n)
y = range(0, stop=1, length=n)
z = range(0, stop=12, length=3)

samplepointssolar = [xi for yj in y, xi in x, zk in z]
samplepointswind = [yj for yj in y, xi in x, zk in z]
samplepointsbat = [zk for yj in y, xi in x, zk in z]

# fill the NaN border cells with the same value

for i = 2:size(curtmatno, 1)
    curtmatno[size(curtmatno,1)-i+2, i] = curtmatno[size(curtmatno,1)-i+1, i]
    curtmat4[size(curtmat4,1)-i+2, i] = curtmat4[size(curtmat4,1)-i+1, i]
    curtmat12[size(curtmat12,1)-i+2, i] = curtmat12[size(curtmat12,1)-i+1, i]
end

curtmat = cat(curtmatno, curtmat4, curtmat12, dims = 3)

# import battery requirements
batteryrequirements = DataFrame(CSV.File("Data/CurtailmentUS/Curtailment_vs_Battery.csv")) |> Matrix
batteryrequirements = batteryrequirements[:, [1, end]]
batteryrequirements = vcat(batteryrequirements, [12 100])
batteryrequirements[:, 2] .= batteryrequirements[:, 2] ./ 100

println("DataLoads.jl has compiled.")
end