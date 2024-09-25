module SteadyStateImp
import Main.SteadyState: p_E_LR
import Main.DataLoads: regionParams

using CSV, DataFrames, LinearAlgebra
import Main.DataLoadsFunc: Matshifter
import Main.Params: majorregions, regions


function solve_steadystate_imp(P, D, M, S, config)
# save original price results
priceresults_orig = copy(p_E_LR)
priceresults_orig_FE = copy(pE_S_FE)

KF_SS = regionParams.KF ./ 10000    # add in minimum of generation
pFK=1
p_F = 10^6

# ---------------------------------------------------------------------------- #
#                               Load New Networks                              #
# ---------------------------------------------------------------------------- #

function create_RegionImp(majorregions, regions)
    mx_of_mx =  [
        Matrix{Float64}(undef, 1474, 726),  # 1474x726 double
        Matrix{Float64}(undef, 1432, 754),  # 1432x754 double
        Matrix{Float64}(undef, 59, 29),     # 59x29 double
        Matrix{Float64}(undef, 114, 52),    # 114x52 double
        Matrix{Float64}(undef, 14, 12),     # 14x12 double
        Matrix{Float64}(undef, 23, 14),     # 23x14 double
        Matrix{Float64}(undef, 76, 45),     # 76x45 double
        Matrix{Float64}(undef, 10, 10),     # 10x10 double
        Matrix{Float64}(undef, 42, 25),     # 42x25 double
        Matrix{Float64}(undef, 157, 77),    # 157x77 double
        Matrix{Float64}(undef, 247, 124),   # 247x124 double
        Matrix{Float64}(undef, 738, 319)    # 738x319 double
    ]

    mutable struct StructRegionImp
        Kmatrix_imp::Matrix{Float64}
        KmatrKmatrix_impix::Matrix{Float64}
        A_imp::Vector{Matrix{Float64}}
        O_imp::Vector{Matrix{Float64}}
        Gam_imp::Vector{Matrix{Float64}}
        Adash_imp::Vector{Matrix{Float64}}
        B_imp::Vector{Matrix{Float64}}
        Zmax_imp::Vector{Float64}

        function StructRegionImp()
            new(
                Matrix{Float64}(undef, 4390, 2532),
                Matrix{Float64}(undef, 4390, 2531),
                mx_of_mx,
                [
                    Matrix{Float64}(undef, 1474, 1474),
                    Matrix{Float64}(undef, 1432, 1432),
                    Matrix{Float64}(undef, 59, 59),
                    Matrix{Float64}(undef, 114, 114),
                    Matrix{Float64}(undef, 14, 14),
                    Matrix{Float64}(undef, 23, 23),
                    Matrix{Float64}(undef, 76, 76),
                    Matrix{Float64}(undef, 10, 10),
                    Matrix{Float64}(undef, 42, 42),
                    Matrix{Float64}(undef, 157, 157),
                    Matrix{Float64}(undef, 247, 247),
                    Matrix{Float64}(undef, 738, 738) 
                ],
                mx_of_mx,
                mx_of_mx,
                [
                    Matrix{Float64}(undef, 726, 726),
                    Matrix{Float64}(undef, 754, 754),
                    Matrix{Float64}(undef, 29, 29),
                    Matrix{Float64}(undef, 52, 52),
                    Matrix{Float64}(undef, 12, 12),
                    Matrix{Float64}(undef, 14, 14),
                    Matrix{Float64}(undef, 45, 45),
                    Matrix{Float64}(undef, 10, 10),
                    Matrix{Float64}(undef, 25, 25),
                    Matrix{Float64}(undef, 77, 77),
                    Matrix{Float64}(undef, 124, 124),
                    Matrix{Float64}(undef, 319, 319)
                ],
                Vector{Float64}(undef, 4390)
            )
        end
    end

    RegionImp = StructRegionImp()



    Kmatrix_imp = DataFrame(CSV.File("Data/ModelDataset/Kmatrix_imp.csv")) |> Matrix
    RegionImp.KmatrKmatrix_impix .= Matshifter(Kmatrix_imp[:, 2:end])

    Linedata_imp = DataFrame(CSV.File("Data/ModelDataset/Linedata_imp.csv"))
    Linecounts_imp = DataFrame(CSV.File("Data/ModelDataset/Linecountsregion_imp.csv"))

    Kmatrix_imp = DataFrame(CSV.File("Data/ModelDataset/Kmatrix_1_imp.csv")) |> Matrix
    Kmatrix_imp = Kmatrix_imp[1:end, 2:majorregions.rowid[1]]

    Omatrix_imp = DataFrame(CSV.File("Data/ModelDataset/Omatrix_1_imp.csv")) |> Matrix
    Omatrix_imp = Omatrix_imp[1:end, 2]

    RegionImp.A_imp[1] .= Matshifter(Kmatrix_imp)
    RegionImp.O_imp[1] .= diagm(Omatrix_imp)
    RegionImp.Gam_imp[1] .= inv(RegionImp.O_imp[1]) * RegionImp.A_imp[1] * inv(RegionImp.A_imp[1]' * inv(RegionImp.O_imp[1]) * RegionImp.A_imp[1])

    indmax = repeat(maximum(RegionImp.Gam_imp[1], dims=2), 1, majorregions.rowid[1] - 1)

    # keep other matrices the same
    for jj = 2:(size(majorregions, 1)-1)
        stringer = "Data/ModelDataset/Kmatrix_$(jj).csv"
        stringer2 = "Data/ModelDataset/Omatrix_$(jj).csv"

        Kmatrix_imp = CSV.read(stringer, DataFrame)
        Kmatrix_imp = Matrix(Kmatrix_imp[1:end, (majorregions.rowid[jj - 1] + 2):majorregions.rowid[jj]])  # MATLAB to Julia index conversion
        Omatrix_imp = CSV.read(stringer2, DataFrame) |> Matrix
        Omatrix_imp = Omatrix_imp[1:end, 2]

        RegionImp.A_imp[jj] .= Matshifter(Kmatrix_imp)
        RegionImp.Adash_imp[jj] .= Matshifter(Kmatrix_imp)
        RegionImp.O_imp[jj] .= diagm(Omatrix_imp)
        RegionImp.Gam_imp[jj] .= inv(RegionImp.O_imp[jj]) * RegionImp.A_imp[jj] * inv(RegionImp.A_imp[jj]' * inv(RegionImp.O_imp[jj]) * RegionImp.A_imp[jj])

    end

    R = size(majorregions, 1) - 1 # regions
    I = 10                        # industries

    Linecounts_imp[!, :rowid2] = [1; Linecounts_imp.rowid[1:end-1].+1]

    # make B matrix
    for jj = 1:size(majorregions, 1)-1
        RegionImp.B_imp[jj] .= inv(RegionImp.A_imp[jj]' * RegionImp.O_imp[jj] * RegionImp.A_imp[jj])
    end

    #max line flows
    RegionImp.Zmax_imp .= Linedata_imp.lcap ./ regions.ffl_capacity_mw[1]

    return RegionImp
end

RegionImp = create_RegionImp(majorregions, regions)

x = range(start = 0.0, stop = 1.0, step = 0.05) 
y = range(start = 0.0, stop = 1.0, step = 0.05) 
z = range(start = 0.0, stop = 12.0, step = 6.0)

imp_curtmat = copy(D.curtmat[:, :, 1:2])

interp2 = LinearInterpolation((x, y), D.curtmat[:,:,1:2])


# ------------------ Solve Power Output within given capital ----------------- #
function solve_power_output_imp(RWParams::StructRWParams, params::StructParams, RunBatteries::Int, RunCurtailment::Int,
                            Initialprod::Int, R_LR::Float64, majorregions::DataFrame, Linecounts::DataFrame, linconscount::Int,
                            regionParams::StructRWParams, curtailmentswitch::Int, interp2, T, mrkteq, config, pB_shifter)
    
    laboralloc_LR, KR_LR_S, KR_LR_W, p_E_LR, w_LR, result_Dout_LR, result_Yout_LR, PC_guess_LR = ss_load_mat()

    global laboralloc_LR
    global KR_LR_S
    global KR_LR_W
    global p_E_LR
    global w_LR
    global result_Dout_LR
    global result_Yout_LR
    global PC_guess_LR
        
    ####
    
    # get guess for long run renewable capital prices
    Deprecations_S = KR_LR_S*params.deltaR
    Deprecations_W = KR_LR_W*params.deltaR
    cumsum_S = sum(Deprecations_S)
    cumsum_W = sum(Deprecations_W)
    cumsum_S = cumsum_S .* (params.iota) .^ (1:500)
    cumsum_W = cumsum_W .* (params.iota) .^ (1:500)
    Qtotal_LR_S = sum(cumsum_S)
    Qtotal_LR_W = sum(cumsum_W)

    Capinvest = ones(1, params.J)
    global Lsector = laboralloc_LR .* params.L
    
    #KR_LR_S = KR_LR_S .+ 0.01
    #KR_LR_W = KR_LR_W .+ 0.01
    KR_LR = KR_LR_S .+ KR_LR_W
    
    if RunBatteries==0 && RunCurtailment==0
        pB_shifter = 0
        config.hoursofstorage=0
    elseif RunCurtailment==1
        pB_shifter = P.pkwh_B / P.pkw_solar
    end  
    
    # output p_B
    p_B = set_battery(KR_LR, config.hoursofstorage, params, Initialprod, T)
    
    
    # initialise run
    #Ddiff = 1
    niters = 1
    niters_in = 1
    diffK = 1.0
    diffp = 1.0
    
    # initiate storage vectors
    sizes = [727, 755, 30, 53, 13, 15, 46, 11, 26, 78, 125, 320, 332]
    last_element = [reshape(Vector{Float64}(undef, sizes[end]), 1, :)]
    
    Lossfac_LR = Matrix{Float64}(undef, 1, 12)
    result_price_LR = [Vector{Float64}(undef, size) for size in sizes[1:end-1]]
    result_price_LR = [result_price_LR..., last_element]
    result_YFout_LR = [Vector{Float64}(undef, size) for size in sizes[1:end-1]]
    result_YFout_LR = [result_YFout_LR..., last_element]
    D_LR = Vector{Float64}(undef, 2531)
    YE_LR = Vector{Float64}(undef, 2531)
    YF_LR = Vector{Float64}(undef, 2531)
    PI_LR = Vector{Float64}(undef, 2531)
    renewshare_LR = zeros(2531, 12)
    SShare_region_LR = zeros(2531, 12)
    curtailmentfactor_S = zeros(1, 2531)
    curtailmentfactor_W = zeros(1, 2531)
    curtailmentfactor = Vector{Float64}(undef, 2531)
    rP_LR = Vector{Float64}(undef, 2531)


    rP_LR .= (R_LR - 1 + params.deltaP).*PC_guess_LR



    while diffK > 10^(-4) 
        println("Number of iterations outer while loop ", niters)
        diffend = 1
        tol = 0.1
        jj=1

        while diffend > tol
            println("Number of iterations inner while loop: ", niters_in)


            # set long run goods prices 
            pg_LR_s = w_LR .^ (params.Vs[:, 1]' .* ones(params.J, 1)) .* 
                    p_E_LR .^ ((params.Vs[:, 2]'.* ones(params.J, 1)) .+ (params.Vs[:, 3]'.* ones(params.J, 1))) .* 
                    ((params.kappa .+ (params.kappa .* p_F_LR ./ p_E_LR) .^ (1 .- params.psi)) .^ (-(params.psi ./ (params.psi .- 1)) .* params.Vs[:, 3]')) .* rP_LR .^ (params.Vs[:, 4]'.* ones(params.J, 1)) ./ 
                    (params.Z .* params.zsector .* params.cdc)
            

            ss_optimize_region_imp!(result_price_LR, result_Dout_LR, result_Yout_LR, result_YFout_LR, Lossfac_LR, pg_LR_s, majorregions, Linecounts, laboralloc_LR, Lsector, params, w_LR,
                rP_LR, p_E_LR, kappa, regionParams, KF_LR, p_F_LR, linconscount, KR_LR_S, KR_LR_W, RegionImp)


            kk = params.N

            KRshifter, YFmax, guess, power, KFshifter, shifter = ss_second_loop(majorregions, Lsector, mrkteq.laboralloc, params, w_LR, rP_LR,
                                                                                result_Dout_LR, result_Yout_LR, pg_LR_s, p_E_LR, kappa,
                                                                                regionParams, KR_LR_S, KR_LR_W, KF_LR, kk, p_F_LR)

            global sub = Matrix{Float64}(undef, majorregions.n[kk], 1)

            Threads.@threads for jj = 1:majorregions.n[kk]
            #for jj = 1:majorregions.n[kk]
                # solve market equilibrium
                local con = [1 -1]
                local guess = [1; KRshifter[jj]]
                local LB = [0; KRshifter[jj]]
                local UB = [10^6; YFmax[jj] + KRshifter[jj]]
                local l_guess = length(guess)

				local x = Vector{Float64}(undef, 2)
                local model = Model(Ipopt.Optimizer)
                set_silent(model)
                @variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
                @constraint(model, c1, con * x <= 0) 
                @objective(model, Min, obj2(x, power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F_LR, params))
                optimize!(model)

                local P_out2 = value.(x)
            
                result_Dout_LR[kk][jj] = P_out2[1]
                result_Yout_LR[kk][jj] = P_out2[2]
                sub[jj] = P_out2[2] - KRshifter[jj]
                result_price_LR[kk][1][1, jj] = Price_Solve(P_out2, shifter[jj], 1, params)[1]
            end

            # used sub to directly define YF_LR below

            for kk=1:(params.N)
                ind = majorregions.rowid2[kk]:majorregions.rowid[kk]
                if kk <= 12
                    p_E_LR[ind] .= result_price_LR[kk]
                    D_LR[ind].= result_Dout_LR[kk]
                    YE_LR[ind] .= result_Yout_LR[kk]
                    YF_LR[ind] .= result_YFout_LR[kk]
                    PI_LR[ind, 1] .= sum(result_price_LR[kk] .* (result_Dout_LR[kk].- result_Yout_LR[kk])) .*
                                        params.L[ind, 1] ./
                                        sum(params.L[ind, 1])
                else
                    result_price_LR[kk][1] .= clamp.(result_price_LR[kk][1], 0.001, 1)
                    matrix_values = result_price_LR[kk][1]
                    p_E_LR[ind] .= vec(matrix_values)
                    D_LR[ind] = result_Dout_LR[kk]
                    YE_LR[ind] = result_Yout_LR[kk]
                    YF_LR[ind] .= sub
                    PI_LR[ind] .= sum(result_price_LR[kk][1] .* (result_Dout_LR[kk].-result_Yout_LR[kk])) .*
                                                                                        params.L[ind, 1] ./
                                                                                        sum(params.L[ind, 1])
                    PI_LR[ind] .= clamp.(PI_LR[ind], 0.0, 0.00001)
                end    
            end

            jj = jj+1

            # get production capital vec
            Ksector = Lsector .* (params.Vs[:,4]'.* ones(params.J, 1)) ./ (params.Vs[:,1]'.* ones(params.J, 1)) .* (w_LR ./ rP_LR)
            global KP_LR = sum(Ksector, dims = 2)

            # update prices and wages
            w_update, w_real, Incomefactor, PC_LR, Xjdashs = wage_update_ms( w_LR, p_E_LR, p_E_LR, p_F_LR, D_LR, YE_LR, rP_LR, KP_LR, PI_LR, 0, params);
            global w_LR = 0.2 .* w_update .+ (1 - 0.2) .* w_LR
            global w_real
            global PC_LR
            # update sectoral allocations
            global laboralloc_LR = Lsector ./ params.L
            relexp = Xjdashs ./ (Xjdashs[:,1] .* ones(1, params.I))
            relab = laboralloc_LR ./ (laboralloc_LR[:,1].* ones(1, params.I))
            global Lsector = Lsector .* clamp.(1 .+ 0.2 .* (relexp .- relab) ./ relab, 0.8, 1.2)
            Lsector = Lsector ./ sum(Lsector, dims=2) .* params.L

            diffend = 0.01
            niters_in += 1
        end

        # update consumption price guess
        PC_guess_LR .= 0.2 .* PC_LR .+ (1 - 0.2) .* PC_guess_LR
        
        # update capital prices
        global Depreciation_LR_S = KR_LR_S .* params.deltaR
        global Depreciation_LR_W = KR_LR_W .* params.deltaR
        Depreciation_LR_S = sum(Depreciation_LR_S)
        Depreciation_LR_W = sum(Depreciation_LR_W)
        cumsum_LR_S = reshape(Depreciation_LR_S .* (params.iota .^ (1:T)), 1, T)
        cumsum_LR_W = reshape(Depreciation_LR_W .* (params.iota .^ (1:T)), 1, T)
        Qtotal_LR_S = sum(cumsum_LR_S)
        Qtotal_LR_W = sum(cumsum_LR_W)
        global p_KR_LR_S = (Initialprod .* (params.iota .^ T) .+ Qtotal_LR_S) .^ (-params.gammaS)
        global p_KR_LR_W = (Initialprod .* (params.iota .^ T) .+ Qtotal_LR_W) .^ (-params.gammaW)

        # get solar shares
        SShare_LR = (regionParams.thetaS ./ p_KR_LR_S) .^ params.varrho ./ ((regionParams.thetaS ./ p_KR_LR_S) .^ params.varrho + (regionParams.thetaW ./ p_KR_LR_W) .^ params.varrho)
        thetabar_LR = regionParams.thetaS .* SShare_LR + regionParams.thetaW .* (1 .- SShare_LR)
        global p_KR_bar_LR = SShare_LR .* p_KR_LR_S + (1 .- SShare_LR) .* p_KR_LR_W

        jj = jj + 1

        # update battery prices
        p_B = update_battery(KR_LR, hoursofstorage, params)
        p_addon = pB_shifter .* hoursofstorage .* p_B

        # generate curtailment factor for grid regions
        for kk in 1:params.N-1
            ind = majorregions.rowid2[kk]:majorregions.rowid[kk]
            global renewshare_LR[1, kk] = 1 - ((sum(YF_LR[ind, :]) ./ sum(YE_LR[ind, :])))
            SShare_region_LR[1, kk] = sum(SShare_LR[ind] .* YE_LR[ind] / sum(YE_LR[ind, :]))
            curtailmentfactor_S[kk] = renewshare_LR[1, kk] * SShare_region_LR[1, kk]
            curtailmentfactor_W[kk] = renewshare_LR[1, kk] * (1 - SShare_region_LR[1, kk])
            
            itp_values = interp2(curtailmentfactor_W[kk], curtailmentfactor_S[kk])
            global curtailmentfactor[ind, 1] .= itp_values
        end

        # generate curtailment factor for off grid-regions
        for kk=params.N
            ind = majorregions.rowid2[kk]:majorregions.rowid[kk]
            global renewshare_LR[ind,1] .= 1 .- YF_LR[ind,:] ./ YE_LR[ind,:]
            SShare_region_LR[ind, 1] = SShare_LR[ind]
            curtailmentfactor_S[ind] .= renewshare_LR[ind, 1] .* SShare_region_LR[ind]
            curtailmentfactor_W[ind] .= renewshare_LR[ind, 1] .* (1 .- SShare_region_LR[ind])
            x_vals = curtailmentfactor_S[ind]
            y_vals = curtailmentfactor_W[ind]
            z_vals = repeat([hoursofstorage], length(ind))
            
            for i in eachindex(ind)
                global curtailmentfactor[ind[i], 1] = interp3(y_vals[i], x_vals[i], z_vals[i])
            end
        end
        
        # if market price is above free entry price, increase renewable capital
        # decrease otherwise

        KR_LR, KR_LR_S, KR_LR_W, diffp, diffK = ss_update_params!(p_KR_bar_LR, p_addon, params, R_LR, regionParams, thetabar_LR, curtailmentswitch,
                                    curtailmentfactor, p_E_LR, KR_LR, KR_LR_S, KR_LR_W, SShare_LR,
                                    diffK, diffp)

        niters += 1
        global KRshifter
    end
    return (KR_LR_S = KR_LR_S, 
            KR_LR_W = KR_LR_W, 
            p_E_LR = p_E_LR, 
            D_LR = D_LR, 
            YE_LR = YE_LR, 
            PC_guess_LR = PC_guess_LR, 
            PI_LR = PI_LR, 
            w_LR = w_LR, 
            laboralloc_LR = laboralloc_LR, 
            p_KR_bar_LR = p_KR_bar_LR, 
            KR_LR = KR_LR, 
            KF_LR = KF_LR, 
            p_KR_LR_S = p_KR_LR_S, 
            p_KR_LR_W = p_KR_LR_W, 
            p_B = p_B,
            p_F_LR = p_F_LR, 
            Lsector = Lsector, 
            YF_LR = YF_LR, 
            diffK = diffK, 
            diffp = diffp, 
            result_Dout_LR = result_Dout_LR, 
            result_Yout_LR = result_Yout_LR, 
            result_YFout_LR = result_YFout_LR,
            result_price_LR = result_price_LR, 
            KP_LR = KP_LR, 
            rP_LR = rP_LR,
            Depreciation_LR_S = Depreciation_LR_S,
            Depreciation_LR_W = Depreciation_LR_W,
            w_real = w_real)
end
sseqI = solve_power_output(RWParams, params, RunBatteries, RunCurtailment,
                                            Initialprod, R_LR, majorregions, Linecounts, linconscount,
                                            regionParams, curtailmentswitch, interp3, T, mrkteq, config, pB_shifter);

println("Steady State diffK= ", sseqI.diffK)
println("Steady State diffp= ", sseqI.diffp)

FUtilization = sseqI.YF_LR ./ sseqI.KF_LR
#mean(FUtilization[1:majorregions.n[1]]))

# Compute electricity and fossil fuel usage in industry and electricity sectors
e2_LR = mrkteq.laboralloc .* repeat(sseqI.D_LR, 1, params.I) .* (repeat(params.Vs[:,2]', params.J, 1) ./ (repeat(params.Vs[:,2]', params.J, 1) + repeat(params.Vs[:,3]', params.J, 1)))
fusage_ind_LR = sum(e2_LR, dims=2) .* (p_E_LR / sseqI.p_F_LR) .^ params.psi
fusage_power_LR = (sseqI.YF_LR ./ sseqI.KF_LR .^ params.alpha2) .^ (1 / params.alpha1)
fusage_total_LR = sum(fusage_power_LR) + sum(fusage_ind_LR)

# compute fossil usage as a share of GDP
GDP = sum(sseqI.w_LR .* params.L .+ p_E_LR .* sseqI.D_LR .+ sseqI.rP_LR .* sseqI.KP_LR .+ sseqI.p_F_LR .* fusage_ind_LR)
Fossshare = sum(sseqI.p_F_LR .* fusage_ind_LR) / GDP

w_real_imp = sseqI.w_LR ./ sseqI.PC_guess_LR
GDP_ind_imp = (sseqI.w_LR .* params.L .+ p_E_LR .* sseqI.D_LR .+ sseqI.rP_LR .* sseqI.KP_LR .+ sseqI.p_F_LR .* fusage_ind_LR) ./ sseqI.PC_guess_LR

# write grid improvement results to Results
priceresults_imp = sseqI.p_E_LR
rmx = [priceresults_imp sseqI.pE_S_FE priceresults_orig priceresults_orig_FE sseq.w_real_LR w_real_imp sseq.GDP_ind_LR GDP_ind_imp params.L]

end
end