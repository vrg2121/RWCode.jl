module MarketFunctions
export calc_subsidyUS, calc_expenditure!, elec_fuel_expenditure!, open_mat_var!,
        market_setup!, second_loop, solve_KP_init!, optimize_region!, optimize_region_test!, solve_initial_equilibrium
using DataFrames, Statistics, MAT
using ..RegionModel, ..MarketEquilibrium
using JuMP, Ipopt

import ..DataLoadsFunc: StructRWParams
import ..ParamsFunctions: StructParams

function solve_initial_equilibrium(params::StructParams, wage_init::Union{Matrix, Vector{Float64}}, majorregions::DataFrame,
    regionParams::StructRWParams, KR_init_S::Matrix{Float64}, KR_init_W::Matrix{Float64}, R_LR::Float64, sectoralempshares::Matrix{Union{Float64, Missing}},
    Linecounts::DataFrame, kappa::Int, regions::DataFrame, linconscount::Int, updw_w::Float64, upw_z::Float64, RWParams::StructRWParams, G::String)
    
    sizes = [727, 755, 30, 53, 13, 15, 46, 11, 26, 78, 125, 320, 332]
    result_price_init = [Vector{Float64}(undef, size) for size in sizes[1:end-1]]
    last_element = [reshape(Vector{Float64}(undef, sizes[end]), 1, :)]
    result_price_init = [result_price_init..., last_element]


    rP_init = Vector{Float64}(undef, 2531)
    pg_init_s = Matrix{Float64}(undef, 2531, 10)
    pE_market_init = Vector{Float64}(undef, 2531)
    Lossfac_init = Matrix{Float64}(undef, 1, 12)
    D_init = Vector{Float64}(undef, 2531)
    YE_init = Vector{Float64}(undef, 2531)
    relexp = Matrix{Float64}(undef, 2531, 10)
    relab = Matrix{Float64}(undef, 2531, 10)
    FUtilization = Vector{Float64}(undef, 2531)	
    Expenditure_init = Vector{Float64}(undef, 2531)
    Lsector = Matrix{Float64}(undef, 2531, 10)

    result_Dout_init = Array{Matrix{Float64}}(undef, 1, 13)
    result_Yout_init = Array{Matrix{Float64}}(undef, 1, 13)
    result_Pout_init = Array{Matrix{Float64}}(undef, 1, 13)

    w_guess, p_E_init, laboralloc, PC_guess_init, p_F_path_guess, 
    wedge, priceshifterupdate, fossilsales = open_mat_var!(result_Dout_init, result_Yout_init, result_Pout_init, G)


    # initial fossil sales guess
    p_F = 0.05

    Lsector .= laboralloc.*params.L

    diffend=1

    power2 = 1 / params.alpha1

        
    while diffend>10^(-2)
        YF_init = Vector{Float64}(undef, 2531)
        PI_init = Vector{Float64}(undef, 2531)
        KP_init = Vector{Float64}(undef, 2531)
        

        # optimization for each region kk
        market_setup!(rP_init, pg_init_s, pE_market_init, wage_init, params, p_E_init, p_F, R_LR, PC_guess_init)

        optimize_region!(result_price_init, result_Dout_init, result_Yout_init, Lossfac_init, majorregions, Linecounts, RWParams, laboralloc,
                            Lsector, params, wage_init, rP_init, linconscount, pg_init_s, pE_market_init, kappa, p_F, regionParams, KR_init_S,
                            KR_init_W);

        # solve market equilibrium
        kk = params.N

        KRshifter, YFmax, guess, power, KFshifter, shifter = second_loop(kk, majorregions, laboralloc, Lsector, params, wage_init, rP_init,
                                                                                result_Pout_init, pg_init_s, pE_market_init, kappa, 
                                                                                regionParams, KR_init_S, KR_init_W, p_F)


        Threads.@threads for jj in 1:majorregions.n[kk]
                local guess = [1; KRshifter[jj]]
                local LB = [0; KRshifter[jj]]
                local UB = [10^6; YFmax[jj] + KRshifter[jj]]
                local l_guess = length(guess)

                local x = Vector{Float64}(undef, l_guess)
                model = Model(Ipopt.Optimizer)
                set_silent(model)
                @variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
                @constraint(model, c1, x[1] - x[2] <= 0) 
                @objective(model, Min, obj2(x, power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F, params))
                optimize!(model)

                local P_out1 = value.(x) 

                if jj == majorregions.n[kk]
                    global P_out = P_out1 
                end 

                local result_Dout_init[kk][jj] = P_out1[1]
                local result_Yout_init[kk][jj] = P_out1[2]

                result_price_init[kk][1][1, jj] = Price_Solve(P_out1, shifter[jj], 1, params)[1]
                
        end

        for kk=1:(params.N)

            ind = majorregions.rowid2[kk]:majorregions.rowid[kk]
            if kk <= 12
                    p_E_init[ind] .= result_price_init[kk]
                    @views D_init[ind].= result_Dout_init[kk]
                    @views YE_init[ind] .= result_Yout_init[kk]
                    @views PI_init[ind, 1] .= (sum(result_price_init[kk] .*(result_Dout_init[kk] - result_Yout_init[kk]))) .*
                                        (params.L[ind, 1]) ./ 
                                        (sum(params.L[ind, 1]))
            
            else
                    result_price_init[kk][1] .= clamp.(result_price_init[kk][1], 0.001, 0.3)	# put bounds on prices
                    matrix_values = result_price_init[kk][1]
                    p_E_init[ind] .= vec(matrix_values)
                    @views D_init[ind].= vec(result_Dout_init[kk])
                    @views YE_init[ind] .= vec(result_Yout_init[kk])
                    @views PI_init[ind, 1] .= (sum(matrix_values.*(result_Dout_init[kk] .- result_Yout_init[kk]))) .*
                                        (params.L[ind, 1]) ./
                                        (sum(params.L[ind, 1]))
            end

        end


        global PI_init

        # get capital vec
        global KP_init = solve_KP_init!(KP_init, Lsector, params, wage_init, rP_init)
        
        # solve wages and calculate productivity
        w_guess
        fossilsales

        w_update, W_Real, Incomefactor, PC_init, Xjdashs = wage_update_ms(w_guess,p_E_init,p_E_init, p_F, D_init, YE_init,rP_init,KP_init,PI_init, fossilsales, params)

        global W_Real
        global PC_init

        w_guess .= 0.1 .* w_update .+ (1 - 0.1) .* w_guess
        params.Z .= params.Z .* clamp.(1 .+ updw_w .* (wage_init .- w_guess) ./ wage_init, 0.95, 1.05)
        diff_w = maximum(abs.(wage_init .- w_guess) ./ wage_init)
        diffend = maximum(diff_w)

        # update relative labor allocations
        laboralloc .= Lsector ./ params.L
        relab .= laboralloc ./ sum(laboralloc, dims=2)
        relexp .= Xjdashs ./ sum(Xjdashs, dims = 2)
        Lsector .= Lsector .* clamp.(1 .+ 0.05 .* (relexp .- relab), 0.9, 1.1)
        Lsector .= Lsector ./ sum(Lsector, dims=2) .* params.L

        # calibrate sectoral z to move sectoral emp closer
        params.zsector[:, :] .= params.zsector[:, :] .* clamp.(1 .+ upw_z .* (sectoralempshares[:, 2:end] .- laboralloc[:, :]), 0.99, 1.01)

        # Update consumption price guess
        PC_guess_init .= 0.2 .* PC_init .+ (1 - 0.2) .* PC_guess_init

        # get fossil fuel usage in the initial steady state
        global @. YF_init = YE_init - regionParams.thetaS * KR_init_S - regionParams.thetaW * KR_init_W
        FUtilization .= YF_init ./ regionParams.maxF

        # calibrate efficiency wedges costs to match initial power prices across countries
        # set subsidy amount for transition
        global subsidy_US = calc_subsidyUS(p_E_init, regions, majorregions) 


        fossilsales, Expenditure_init = elec_fuel_expenditure!(laboralloc, D_init, params, p_E_init, p_F, YF_init, regionParams,
                                                                        wage_init, rP_init, KP_init, fossilsales, YE_init, regions, PI_init)

        println("diffend = ", diffend)
    end

    return (Lsector = Lsector, 
            result_price_init = result_price_init, 
            pg_init_s = pg_init_s, 
            pE_market_init = pE_market_init, 
            Lossfac_init = Lossfac_init, 
            D_init = D_init, 
            YE_init = YE_init, 
            diffend = diffend, 
            w_guess = w_guess, 
            p_E_init = p_E_init, 
            result_Pout_init = result_Pout_init, 
            laboralloc = laboralloc, 
            PC_guess_init = PC_guess_init, 
            result_Dout_init = result_Dout_init, 
            result_Yout_init = result_Yout_init, 
            p_F_path_guess = p_F_path_guess, 
            wedge = wedge, 
            priceshifterupdate = priceshifterupdate, 
            fossilsales = fossilsales, 
            P_out = P_out, 
            Expenditure_init = Expenditure_init,
            rP_init = rP_init, 
            PI_init = PI_init, 
            KP_init = KP_init, 
            W_Real = W_Real, 
            PC_init = PC_init, 
            YF_init = YF_init, 
            p_F = p_F,
            subsidy_US = subsidy_US)
end

# ---------------------------------------------------------------------------- #
#                      Functions to Solve the Equilibrium                      #
# ---------------------------------------------------------------------------- #

function open_mat_var!(result_Dout_init::Matrix{Matrix{Float64}}, result_Yout_init::Matrix{Matrix{Float64}}, result_Pout_init::Matrix{Matrix{Float64}}, G::String)

    wf = matopen("$G/w_guess_mat.mat")
    w_guess = read(wf, "w_guess")
    close(wf)

    pEf = matopen("$G/p_E_guessmat.mat") 
    p_E_init = read(pEf, "p_E_init")
    close(pEf)
    p_E_init = vec(p_E_init)

    pof = matopen("$G/Pout_guess_init.mat")
    result_Pout_init .= read(pof, "result_Pout_init")
    close(pof)

    lf = matopen("$G/laboralloc_guess.mat")
    laboralloc = read(lf, "laboralloc")
    close(lf)

    PCf = matopen("$G/PC_guess_init.mat")
    PC_guess_init = read(PCf, "PC_guess_init")
    close(PCf)

    dof = matopen("$G/Dout_guess_init.mat")
    result_Dout_init .= read(dof, "result_Dout_init")
    close(dof)

    yof = matopen("$G/Yout_guess_init.mat")
    result_Yout_init .= read(yof, "result_Yout_init")
    close(yof)

    pff = matopen("$G/p_F_path_guess_saved.mat")
    p_F_path_guess = read(pff, "p_F_path_guess")
    close(pff)

    wf = matopen("$G/wedge_vec.mat")
    wedge = read(wf, "wedge")
    close(wf)

    psf = matopen("$G/priceshifterupdate_vec.mat")
    priceshifterupdate = read(psf, "priceshifterupdate")
    close(psf)

    ff = matopen("$G/fossilsales_guess.mat")
    fossilsales = read(ff, "fossilsales")
    close(ff)

    return w_guess, p_E_init, laboralloc, PC_guess_init, p_F_path_guess, wedge, priceshifterupdate, fossilsales
end

function market_setup!(rP_init::Vector, pg_init_s::Matrix, pE_market_init::Vector, 
    wage_init::Vector, params::StructParams, p_E_init::Vector, p_F::Float64, R_LR::Float64, PC_guess_init::Matrix)
    
    # Set initial capital returns
    rP_init .= (R_LR - 1 + params.deltaP) .* PC_guess_init      # long run return on production capital

    # Calculate final good price
    pg_init_s .= (wage_init .^ (params.Vs[:, 1]' .* ones(params.J, 1))) .* 
                (p_E_init .^ ((params.Vs[:, 2]' .* ones(params.J, 1)) .+ (params.Vs[:, 3]' .* ones(params.J, 1)))) .* 
                ((params.kappa .+ (params.kappa .* p_F ./ p_E_init) .^ (1 - params.psi)) .^ 
                (-(params.psi / (params.psi - 1)) .* params.Vs[:, 3]')) .* 
                (rP_init .^ (params.Vs[:, 4]' .* ones(params.J, 1))) ./ 
                (params.Z .* params.zsector .* params.cdc)

    pE_market_init .= copy(p_E_init)

end


function optimize_region!(result_price_init::Vector, result_Dout_init::Matrix{Matrix{Float64}}, result_Yout_init::Matrix{Matrix{Float64}}, Lossfac_init::Matrix,
	majorregions::DataFrame, Linecounts::DataFrame, RWParams::StructRWParams, laboralloc::Matrix, Lsector::Matrix, params::StructParams,
	wage_init::Vector, rP_init::Vector, linconscount::Int, pg_init_s::Matrix, pE_market_init::Vector, kappa::Int,
	p_F::Float64, regionParams::StructRWParams, KR_init_S::Matrix, KR_init_W::Matrix)

    # set up optimization problem for region kk
	Threads.@threads :static for kk in 1:(params.N - 1)
		
        # the local call in this section localizes the memory to each thread to reduce crossing data 
		local l_guess, LB, UB, guess, power, shifter, KFshifter, KRshifter, n, mid = data_set_up(kk, majorregions, Linecounts, RWParams, laboralloc,
											Lsector, params, wage_init, rP_init, pg_init_s, pE_market_init, kappa, regionParams, 
											regionParams.KF, p_F, linconscount, KR_init_S, KR_init_W, "market")

		# solve the model for region kk (see RegionModel.jl)
		local P_out = solve_model(kk, l_guess, LB, UB, guess, regionParams, params, power, shifter, KFshifter, KRshifter, p_F, mid)

		result_price_init[kk].=Price_Solve(P_out, shifter, n, params) #.MarketEquilibrium.jl
		@views result_Dout_init[kk] .= P_out[1:mid]
		@views result_Yout_init[kk] .= P_out[1+mid:end]
		local Pvec = P_out[1+mid:end] .- P_out[1:mid]
		local Losses = Pvec[2:end]' * regionParams.B[kk] * Pvec[2:end] .* params.Rweight
		@views Lossfac_init[kk] = Losses ./ sum(result_Yout_init[kk])

	end
end


function second_loop(kk::Int, majorregions::DataFrame, laboralloc::Matrix, Lsector::Matrix, params::StructParams, 
            wage_init::Vector, rP_init::Vector, result_Pout_init::Matrix, 
            pg_init_s::Matrix, pE_market_init::Vector, kappa::Int, 
            regionParams::StructRWParams, KR_init_S::Matrix, KR_init_W::Matrix, p_F::Float64)

    YFmax = Vector{Float64}(undef, 332)
    KRshifter = Vector{Float64}(undef, 332)
    guess = Matrix{Float64}(undef, 348, 1)
    power = Matrix{Float64}(undef, 332, 10)
    KFshifter = Vector{Float64}(undef, 332)
    shifter = Matrix{Float64}(undef, 332, 10)
    prices = Matrix{Float64}(undef, 332, 10)
    pg_s = Matrix{Float64}(undef, 332, 10)
    Lshifter = Matrix{Float64}(undef, 332, 10)
    secalloc = Matrix{Float64}(undef, 332, 10)
    Kshifter = Matrix{Float64}(undef, 332, 10)


    # set up optimization problem for region kk
    ind = majorregions.rowid2[kk]:majorregions.rowid[kk]

    @views secalloc .= laboralloc[ind, :]
    @views Lshifter .= Lsector[ind, :]
    @views Kshifter .= Lshifter .* 
        (params.Vs[:, 4]' .* ones(majorregions.n[kk], 1)) ./ 
        (params.Vs[:, 1]' .* ones(majorregions.n[kk], 1)) .* 
        wage_init[ind] ./ 
        rP_init[ind]

    Jlength::Int = majorregions.n[kk]
    guess .= result_Pout_init[kk]

    # define shifters for objective function
    @views pg_s .= pg_init_s[ind,:]
    @views prices .= (pE_market_init[ind] .* ones(1, params.I))
    @views power .= (params.Vs[:, 2]' .* ones(Jlength, 1))+(params.Vs[:,3]' .* ones(Jlength, 1))
    shifter .= pg_s .* 
            ((kappa .+ (prices ./ (kappa .* p_F))) .^ (params.psi - 1)) .^ 
            ((params.psi / (params.psi - 1)) .* (params.Vs[:, 3]' .* ones(majorregions.n[kk], 1))) .* 
            (1 .+ (params.Vs[:, 3]' .* ones(majorregions.n[kk], 1)) ./ (params.Vs[:, 2]' .* ones(majorregions.n[kk], 1))) .^ 
            (-(params.Vs[:, 2]' .* ones(majorregions.n[kk], 1)) .- (params.Vs[:, 2]' .* ones(majorregions.n[kk], 1))) .* 
            params.Z[ind] .* params.zsector[ind, :] .* 
            (Lshifter .^ (params.Vs[:, 1]' .* ones(majorregions.n[kk], 1))) .* 
            (Kshifter .^ (params.Vs[:, 4]' .* ones(majorregions.n[kk], 1)))

    shifter .= shifter .* secalloc .^ power

    @views KRshifter .= regionParams.thetaS[ind] .* KR_init_S[ind] +
                        regionParams.thetaW[ind] .* KR_init_W[ind]
    @views KFshifter .= regionParams.KF[ind]

    # define bounds
    @views YFmax .= regionParams.maxF[ind]

    # define shifter
    return KRshifter, YFmax, guess, power, KFshifter, shifter

end

function solve_KP_init!(KP_init::Vector, Lsector::Matrix, params::StructParams, wage_init::Vector, rP_init::Vector)
    # get capital vec
    Ksector = Lsector .* 
            (params.Vs[:, 4]' .* ones(params.J, 1)) ./
            (params.Vs[:, 1]' .* ones(params.J, 1)) .* 
            (wage_init ./ rP_init)
    KP_init .= sum(Ksector, dims = 2)
    return KP_init
end

function calc_subsidyUS(p_E_init::Vector, regions::DataFrame, majorregions::DataFrame)
    p_E_weighted = Vector{Float64}(undef, 2531)
    @. p_E_weighted = p_E_init * regions.employment / regions.employment_region
    pricetable = DataFrame(p_E_weighted = p_E_weighted, region_groups = regions.region_groups)
    tblstats = combine(groupby(pricetable, :region_groups), :p_E_weighted => mean => :mean_p_E_weighted)
    p_E_weighted = tblstats.mean_p_E_weighted .* majorregions.n

    # set subsidy amount for transition
    subsidy_US=p_E_weighted[1]*0.15
    return subsidy_US
end

function elec_fuel_expenditure!(laboralloc::Matrix, D_init::Vector, params::StructParams, p_E_init::Vector, p_F::Float64, 
    YF_init::Vector, regionParams::StructRWParams, wage_init::Vector, rP_init::Vector, KP_init::Vector, fossilsales::Matrix, 
    YE_init::Vector, regions::DataFrame, PI_init::Vector)
    e2_init = Matrix{Float64}(undef, 2531, 10)
    fusage_ind_init = Vector{Float64}(undef, 2531)
    fusage_power_init = Vector{Float64}(undef, 2531)
    Expenditure_init = Vector{Float64}(undef, 2531)

    calc_expenditure!(e2_init, fusage_ind_init, fusage_power_init, Expenditure_init, fossilsales, 
    laboralloc, D_init, params, p_E_init, p_F, YF_init, regionParams,
    wage_init, rP_init, KP_init, YE_init, regions, PI_init)

    return fossilsales, Expenditure_init
end

function calc_expenditure!(e2_init::Matrix, fusage_ind_init::Vector, 
    fusage_power_init::Vector, Expenditure_init::Vector,
    fossilsales::Matrix, laboralloc::Matrix, D_init::Vector, params::StructParams, p_E_init::Vector, p_F::Float64, 
    YF_init::Vector, regionParams::StructRWParams, wage_init::Vector, rP_init::Vector, KP_init::Vector, YE_init::Vector,
    regions::DataFrame, PI_init::Vector)
    
    # compute electricty and fossil fuel usage in industry and electricity
    e2_init .= laboralloc .* (D_init .* ones(1, params.I)) .* ((params.Vs[:,2]' .* ones(params.J, 1)) ./ 
    ((params.Vs[:,2]' .* ones(params.J, 1))+(params.Vs[:,3]' .* ones(params.J, 1))))
    fusage_ind_init .= (params.kappa .^ -1) .* sum(e2_init, dims=2) .* (p_E_init ./ p_F) .^ params.psi
    fusage_power_init .= (YF_init ./ regionParams.maxF .^ params.alpha2) .^ (1 / params.alpha1)
    fusage_total_init = sum(fusage_power_init) + sum(fusage_ind_init)

    # compute fossil usage as a share of GDP
    #GDP=sum(wage_init.*params.L+p_E_init.*D_init+rP_init.*KP_init+p_F.*fusage_ind_init)
    #Fossshare_init=sum(p_F.*fusage_total_init)./GDP

    # generate expenditure shares on each
    fossilsales .= p_F .* fusage_total_init .* regions.reserves_share
    Expenditure_init .= wage_init .* params.L + rP_init .* KP_init + fossilsales + 
                        p_E_init .* YE_init - p_F .* fusage_power_init + p_F .* (fusage_power_init + fusage_ind_init) + PI_init
    #fossshare_j = fossilsales ./ Expenditure_init
end


end