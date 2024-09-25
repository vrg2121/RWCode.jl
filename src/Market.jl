# Market.jl
module Market

# export variables
export solve_market

# load functions
using Main.DataAdjustments
using Main.MarketFunctions
using Main.MarketEquilibrium
using Main.RegionModel
import Main.ModelConfiguration: ModelConfig

# load packages
using JuMP, Ipopt
import Random: Random
import LinearAlgebra: I
import SparseArrays: sparse
import JLD2: @save
import DataFrames: DataFrame

function solve_market(P::NamedTuple, DL::NamedTuple, config::ModelConfig, G::String)
# ---------------------------------------------------------------------------- #
#                           Solve Market Equilibrium                           #
# ---------------------------------------------------------------------------- #


    mrkteq = solve_initial_equilibrium(P.params, DL.wage_init, P.majorregions,
                                        DL.regionParams, DL.KR_init_S, DL.KR_init_W, DL.R_LR, DL.sectoralempshares,
                                        P.Linecounts, P.kappa, P.regions, P.linconscount, P.updw_w, P.upw_z, DL.RWParams, G);
    println("Initial Calibration= ", mrkteq.diffend)

    Z = P.params.Z
    Zsec = P.params.zsector
    w_guess = mrkteq.w_guess
    p_E_init = mrkteq.p_E_init
    result_Dout_init = mrkteq.result_Dout_init
    result_Yout_init = mrkteq.result_Yout_init
    PC_guess_init = mrkteq.PC_guess_init
    laboralloc = mrkteq.laboralloc
    wedge = mrkteq.wedge
    priceshifterupdate = mrkteq.priceshifterupdate
    fossilsales = mrkteq.fossilsales
    subsidy_US = mrkteq.subsidy_US


    # Save variables
    @save "$G/w_guess_mat.jld2" w_guess
    @save "$G/p_E_guessmat.jld2" p_E_init
    @save "$G/Dout_guess_init.jld2" result_Dout_init
    @save "$G/Yout_guess_init.jld2" result_Yout_init
    @save "$G/PC_guess_init.jld2" PC_guess_init

    @save "$G/laboralloc_guess.jld2" laboralloc
    @save "$G/z_mat.jld2" Z
    @save "$G/z_sec_mat.jld2" Zsec
    @save "$G/wedge_vec.jld2" wedge
    @save "$G/priceshifterupdate_vec.jld2" priceshifterupdate
    @save "$G/fossilsales_guess.jld2" fossilsales

    # set initial power output vector
    P_out_init = mrkteq.P_out

    # set Q initial to be the solar already installed
    Qtotal_init_S = sum(DL.KR_init_S)
    Qtotal_init_W = sum(DL.KR_init_W)
    p_KR_init_S=(config.Initialprod+Qtotal_init_S).^(-P.params.gammaS);    
    p_KR_init_W=(config.Initialprod+Qtotal_init_W).^(-P.params.gammaW);
    SShare_init=(DL.regionParams.thetaS./p_KR_init_S).^P.params.varrho./((DL.regionParams.thetaS./p_KR_init_S).^P.params.varrho+(DL.regionParams.thetaW./p_KR_init_W).^P.params.varrho)
    thetabar_init = DL.regionParams.thetaS .* SShare_init + DL.regionParams.thetaW .* (1 .- SShare_init)
    p_KR_bar_init = SShare_init .* p_KR_init_S + (1 .- SShare_init) .* p_KR_init_W
    pE_FE_init=(p_KR_bar_init-p_KR_bar_init*(1-P.params.deltaR)./DL.R_LR).*DL.regionParams.costshifter./thetabar_init

    wageresults = Matrix{Float64}(undef, 2531, 2)
    priceresults = Vector{Float64}(undef, 2531)
    PCresults = Vector{Float64}(undef, 2531)
    KF_init = Vector{Float64}(undef, 2531)

    wageresults[:,1].=copy(mrkteq.W_Real)
    priceresults[:,1].=copy(p_E_init)
    laboralloc_init=copy(laboralloc)
    PCresults[:,1].=copy(p_E_init)


    #get renewable shares
    renewshareUS=1-(sum(mrkteq.YF_init[1:P.majorregions.n[1],:])./sum(mrkteq.YE_init[1:P.majorregions.n[1]]));
    renewshareEU=1-(sum(mrkteq.YF_init[P.majorregions.rowid[1]+1:P.majorregions.rowid[2]]))./sum(mrkteq.YE_init[P.majorregions.rowid[1]+1:P.majorregions.rowid[2]])

    p_F_int=copy(mrkteq.p_F)
    KF_init=copy(DL.regionParams.KF)

    return (
        wageresults = wageresults,
        p_KR_bar_init = p_KR_bar_init,
        KF_init = KF_init,
        laboralloc_init = laboralloc_init,
        p_KR_init_S = p_KR_init_S,
        p_KR_init_W = p_KR_init_W,
        renewshareUS = renewshareUS,
        p_F_int = p_F_int,
        mrkteq = mrkteq,
        priceresults = priceresults
)

end
end