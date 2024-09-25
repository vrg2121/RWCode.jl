module SteadyState

using Main.Params, MAT, Main.ParameterizeRun, Main.DataLoads, Main.DataAdjustments, Main.MarketEquilibrium, Main.Market, Printf, Random, LinearAlgebra, SparseArrays, JuMP, Ipopt, JLD2, Statistics, DataFrames, JLD2, Interpolations, Plots, SparseArrays
using Main.DataLoads: samplepointssolar, samplepointsbat, samplepointswind
using Main.SteadyStateFunction
export KR_LR_S, KR_LR_W, p_E_LR, D_LR, YE_LR, PC_guess_LR, PI_LR, w_LR, laboralloc_LR, p_KR_bar_LR, KR_LR, KF_LR, p_KR_LR_S, p_KR_LR_W, p_B,
       p_F_LR, Lsector, YF_LR, curtailmentfactor, interp3, GDP, wr, wagechange, welfare_wagechange, welfare_capitalchange, 
       welfare_electricitychange, welfare_fossilchange, KRshfiter, sub, P_out


function isint(x)
    x == round(x)
end

# ---------------------------------------------------------------------------- #
#                             Step 0: Initial Guess                            #
# ---------------------------------------------------------------------------- #


"""
NOTES ON THE TRANSLATION FROM JULIA TO MATLAB

Manually Solving NLP Optimization Problem using Ipopt.jl
user defined gradient:
- Ipopt optimizer is unstable when it comes to multithreading and throws errors with large numbers of threads
- the first multithreaded loop keeps ending up with a negative number in x, which is technically impossible
    - checking the error suggests that there is instability in the model that is difficult to resolve
    - the error only crops up after the while loops have run 6-22 times
- A user defined gradient slows down the model because you cannot multithread it without creating errors
- results from model with user defined gradient the same as results from automatic solve

user defined hessian:
- There is not a good way to extract the lagrangian as the model is being solved in order to manually calculate the hessian



figuring out how to calculate the Interpolation
    MATLAB:

    curtailmentfactor(majorregions.rowid2(kk):majorregions.rowid(kk),1)=interp3(samplepointssolar,samplepointswind,samplepointsbat, curtmat,curtailmentfactor_S(kk), curtailmentfactor_W(kk), hoursofstorage);
    
    
    linear interpolation
    samplepointssolar, samplepointswind, samplepointsbat = coordinates of the sample points
    curtmat = corresponding function values at each sample points
    curtailmentfactor_S(kk), curtailmentfactor_W(kk), hoursofstorage = coordinates of query points

    because the samplepoints are initiated using meshgrid, this is actually a scattered data interpolation

TODO: 
- reset local variables for the multithreading
"""
Expenditure_init

laboralloc_LR, KR_LR_S, KR_LR_W, p_E_LR, w_LR, result_Dout_LR, 
                                    result_Yout_LR, PC_guess_LR = ss_load_mat()


# ------------------ Solve Power Output within given capital ----------------- #
function solve_power_output(RWParams::Main.DataLoads.StructRWParams, params::Main.Params.StructParams, RunBatteries, RunCurtailment,
                            Initialprod, R_LR, majorregions::DataFrame, Linecounts::DataFrame, linconscount,
                            regionParams::Main.DataLoads.StructRWParams, curtailmentswitch)    

    global laboralloc_LR
    global KR_LR_S
    global KR_LR_W
    global p_E_LR
    global w_LR
    global result_Dout_LR
    global result_Yout_LR
    global PC_guess_LR

    KF_LR=RWParams.KF/10000    # add in minimum of generation
    #maxF_SS=KF_LR
    #pFK=1                  # price of fossil fuel capital
    
    p_F_LR=1
    
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
    
    global Lsector = laboralloc_LR .* params.L
    
    global KR_LR_S = KR_LR_S .+ 0.01
    global KR_LR_W = KR_LR_W .+ 0.01
    global KR_LR = KR_LR_S .+ KR_LR_W
    
    if RunBatteries==0 && RunCurtailment==0
        pB_shifter = 0
        hoursofstorage=0
    end  
    
    # output p_B
    p_B = set_battery(KR_LR, hoursofstorage, params, Initialprod, T)
    
    
    # initialise run
    #Ddiff = 1
    niters = 1
    niters_in = 1
    global diffK = 1.0
    global diffp = 1.0
    
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

    #using LazyGrids

    x = range(start = 0.0, stop = 1.0, step = 0.05) 
    y = range(start = 0.0, stop = 1.0, step = 0.05) 
    z = range(start = 0.0, stop = 12.0, step = 6.0)

    interp = interpolate((x, y, z), curtmat, Gridded(Linear()))
    
    #x = unique(samplepointssolar[:, :, 1])
    #y = unique(samplepointswind[:, :, 1])   
    #z = unique(samplepointsbat[1, 1, :]) 
    #interp3 = interpolate((samplepointssolar, samplepointswind, samplepointsbat), curtmat, Gridded(Linear()))


    global rP_LR = (R_LR - 1 + params.deltaP)*PC_guess_LR



    while diffK > 10^(-2)
     
        println("Number of iterations outer while loop: ", niters)
        diffend = 1
        #tol = 0.01
        jj=1

        #while diffend > tol
            #println("Number of iterations inner while loop: ", niters_in)


            # set long run goods prices 
            pg_LR_s = w_LR .^ (params.Vs[:, 1]' .* ones(params.J, 1)) .* 
                    p_E_LR .^ ((params.Vs[:, 2]'.* ones(params.J, 1)) .+ (params.Vs[:, 3]'.* ones(params.J, 1))) .* 
                    ((params.kappa .+ (params.kappa .* p_F_LR ./ p_E_LR) .^ (1 .- params.psi)) .^ (-(params.psi ./ (params.psi .- 1)) .* params.Vs[:, 3]')) .* rP_LR .^ (params.Vs[:, 4]'.* ones(params.J, 1)) ./ 
                    (params.Z .* params.zsector .* params.cdc)

            Threads.@threads for kk in 1:(params.N - 1)
            #for kk = 1:(params.N - 1)
                local ind = majorregions.rowid2[kk]:majorregions.rowid[kk]
                local n = majorregions.n[kk]
			    local l_ind = Linecounts.rowid2[kk]:Linecounts.rowid[kk]
                local gam = RWParams.Gam[kk]
			    local l_n = Linecounts.n[kk]

                # set up optimization problem for region kk
                local @views secalloc = laboralloc_LR[ind, :]
                local @views Lshifter = Lsector[ind, :]
                local Kshifter = Lsector[ind, :] .* 
                        (params.Vs[:, 4]'.* ones(n, 1)) ./(params.Vs[:, 1]'.* ones(n, 1)) .* 
                        w_LR[ind] ./rP_LR[ind]
                #Ltotal = sum(Lshifter, dims=2)
                local Jlength = n

                # define data for inequality constraints
                local linecons = RWParams.Zmax[l_ind]
                local Gammatrix = [zeros(size(gam, 1)) gam]
                if linconscount < l_n
                    Random.seed!(1)
                    randvec = rand(l_n)
                    randvec = randvec .> (linconscount / l_n)
                    Gammatrix[randvec, :] .= 0
                end
                local stacker = [-Matrix(I, n, n) Matrix(I, n, n)]
                local Gammatrix = sparse(Gammatrix * stacker)
                local Gammatrix = [Gammatrix; -Gammatrix]
                local linecons = [linecons; linecons]

                # define shifters for objective function
                local pg_s = pg_LR_s[ind, :]
                local prices = (p_E_LR[ind] .* ones(1, params.I))
                local power = ((params.Vs[:, 2]'.* ones(Jlength, 1)) .+ (params.Vs[:, 3]'.* ones(Jlength, 1)))
                
                local shifter = pg_s .* (kappa .+ (prices ./ (kappa .* p_F_LR)).^(params.psi - 1)) .^ ((params.psi / (params.psi - 1)) .* (params.Vs[:, 3]'.* ones(n, 1))) .* 
                          (1 .+ (params.Vs[:,3]'.* ones(n, 1)) ./ (params.Vs[:,2]'.* ones(n, 1))) .^ (-(params.Vs[:,2]'.* ones(n, 1)) - 
                          (params.Vs[:,2]'.* ones(n, 1))) .* 
                          params.Z[ind] .* 
                          params.zsector[ind, :] .*
                          Lshifter .^ (params.Vs[:,1]'.* ones(n, 1)) .* 
                          Kshifter .^ (params.Vs[:,4]'.* ones(n, 1))

                local shifter=shifter.*secalloc.^power
                local KRshifter = regionParams.thetaS[ind] .* KR_LR_S[ind] +
                            regionParams.thetaW[ind] .* KR_LR_W[ind]
                local KFshifter = KF_LR[ind]

                # define bounds
                local YFmax =KF_LR[ind]
                local LB = [zeros(n); KRshifter]
                if any(LB .< 0)
                    println("Value in LB is negative, iteration $niters, thread $kk")
                end
                local UB = [1e3 * ones(n); YFmax .+ KRshifter .+ 1]

                # define bounds
                local guess = [KRshifter; KRshifter .+ 0.001]
                local l_guess = length(guess)

                #solve market equilibrium
                local x = Vector{Float64}(undef, l_guess)
                local model = Model(Ipopt.Optimizer);
                set_silent(model)
                #set_attribute(model, "max_iter", 10^5)
                @variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
                @constraint(model, c1, mycon(x, regionParams.B[kk], params) == 0)

                #function new_obj_f(x...)
                #    return obj(collect(x), power, shifter, KFshifter, KRshifter, p_F_LR, params)
                #end
                
                #function new_grad_f(g, x...)
                #    g = grad_f(collect(x), power, shifter, KFshifter, KRshifter, p_F_LR, params)
                #    return
                #end

                #function new_hess(h, x...)
                #    h = hessinterior(collect(x), l, power, shifter, KFshifter, KRshifter, p_F_LR, regionParams.B[kk], params)
                #    return
                #end

                @objective(model, Min, obj(x, power, shifter, KFshifter, KRshifter, p_F_LR, params))
                #@operator(model, mrkt_eq, l_guess, new_obj_f, new_grad_f, new_hess)
                #@objective(model, Min, mrkt_eq(x...))
                optimize!(model)

                #l = dual(c1)

                local P_out = value.(x)
                #global P_out1 = P_out
                result_price_LR[kk] .= Price_Solve(P_out, shifter, Jlength, params)
                result_Dout_LR[kk] .= P_out[1:end÷2]
                result_Yout_LR[kk] .= P_out[end÷2+1:end]
                result_YFout_LR[kk] .= P_out[end÷2+1:end] .- KRshifter
                local Pvec = P_out[end÷2+1:end] .- P_out[1:end÷2]
                local Losses = Pvec[2:end]' * regionParams.B[kk] * Pvec[2:end] .* params.Rweight
                Lossfac_LR[1, kk] = Losses / sum(result_Yout_init[kk])
            end

            kk = params.N

            KRshifter, YFmax, guess, power, KFshifter, shifter = ss_second_loop(majorregions, Lsector, laboralloc, params, w_LR, rP_LR,
                                                                                result_Dout_LR, result_Yout_LR, pg_LR_s, p_E_LR, kappa,
                                                                                regionParams, KR_LR_S, KR_LR_W, KF_LR, kk, p_F_LR)

            sub = Matrix{Float64}(undef, majorregions.n[kk], 1)

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
                #function new_obj2(x...)
                #    return obj2(collect(x), power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F_LR, params)
                #end
                
                #function new_grad2(g, x...)
                #    g = grad_f(collect(x), power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F_LR, params)
                #    return
                #end
                #@operator(model, mrkt_eq2, l_guess, new_obj2, new_grad2)
                #@objective(model, Min, mrkt_eq2(x...))
                @objective(model, Min, obj2(x, power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F_LR, params))
                optimize!(model)

                local P_out2 = value.(x)

                if jj == majorregions.n[kk]
					global P_out = P_out2  # Only set to global on specific condition
				end 
                
                result_Dout_LR[kk][jj] = P_out2[1]
                result_Yout_LR[kk][jj] = P_out2[2]
                
                sub[jj] = P_out2[2] - KRshifter[jj]

                #local result_price = Price_Solve(P_out2, shifter[jj], 1, params)[1]
                #local result_price = clamp(result_price, 0.001, 1)
                #result_price_LR[kk][1][1, jj]  = result_price
                
                result_price_LR[kk][1][1, jj] = Price_Solve(P_out2, shifter[jj], 1, params)[1]
                #result_price_LR[kk][1][1, jj] = clamp(result_price_LR[kk][1][1, jj], 0.001, 1)
            end

            #global result_price_LR[13][1] .= clamp.(result_price_LR[13][1], 0.001, 1)
            
            #for jj = 1:majorregions.n[kk]
            #    result_price_LR[kk][1][1, jj] = max(result_price_LR[kk][1][1, jj], 0.001)
            #    result_price_LR[kk][1][1, jj] = min(result_price_LR[kk][1][1, jj], 1)
            #end

            for jj = 1:majorregions.n[kk]
                result_YFout_LR[kk][1][jj] = sub[jj, 1]
            end

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
                    YF_LR[ind] = result_YFout_LR[kk][1]
                    PI_LR[ind] .= sum(result_price_LR[kk][1] .* (result_Dout_LR[kk].-result_Yout_LR[kk])) .*
                                                                                        params.L[ind, 1] ./
                                                                                        sum(params.L[ind, 1])
                end    
            end

            jj = jj+1

            # get production capital vec
            Ksector = Lsector .* (params.Vs[:,4]'.* ones(params.J, 1)) ./ (params.Vs[:,1]'.* ones(params.J, 1)) .* (w_LR ./ rP_LR)
            global KP_LR = sum(Ksector, dims = 2)

            # update prices and wages
            w_update, w_real, Incomefactor, PC_LR, Xjdashs = wage_update_ms( w_LR, p_E_LR, p_E_LR, p_F_LR, D_LR, YE_LR, rP_LR, KP_LR, PI_LR, 0, params);
            #diff_w = maximum(abs.(w_LR .- w_update) ./ w_LR)
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
        global PC_guess_LR .= 0.2 .* PC_LR .+ (1 - 0.2) .* PC_guess_LR
        
        # update capital prices
        global Depreciation_LR_S = KR_LR_S .* params.deltaR
        global Depreciation_LR_W = KR_LR_W .* params.deltaR
        global Depreciation_LR_S = sum(Depreciation_LR_S)
        global Depreciation_LR_W = sum(Depreciation_LR_W)
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
            
            #itp_values = [interp3(curtailmentfactor_S[kk], curtailmentfactor_W[kk], hoursofstorage) for n in ind]
            itp_values = interp(curtailmentfactor_S[kk], curtailmentfactor_W[kk], hoursofstorage)
            global curtailmentfactor[ind, 1] .= itp_values
                                             
        #end

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
                #global curtailmentfactor[ind[i], 1] = interp3[x_vals[i], y_vals[i], z_vals[i]] 
                global curtailmentfactor[ind[i], 1] = interp(x_vals[i], y_vals[i], z_vals[i])
            end
        end
        
        # if market price is above free entry price, increase renewable capital
        # decrease otherwise

        global KR_LR, KR_LR_S, KR_LR_W, diffp, diffK = ss_update_params!(p_KR_bar_LR, p_addon, params, R_LR, regionParams, thetabar_LR, curtailmentswitch,
                                    curtailmentfactor, p_E_LR, KR_LR, KR_LR_S, KR_LR_W, SShare_LR,
                                    diffK, diffp)

        println("Steady State diffK= $diffK")
        println("Steady State diffp= $diffp")
                                                        
        global KRshifter

        niters += 1
    end
    return KR_LR_S, KR_LR_W, p_E_LR, D_LR, YE_LR, PC_guess_LR, PI_LR, w_LR, laboralloc_LR, p_KR_bar_LR, KR_LR, KF_LR, p_KR_LR_S, p_KR_LR_W, p_B,
        p_F_LR, Lsector, YF_LR, curtailmentfactor, diffK, diffp, KRshifter, P_out, result_Dout_LR, result_Yout_LR, result_YFout_LR,
        result_price_LR, KP_LR
end
    KR_LR_S, KR_LR_W, p_E_LR, D_LR, YE_LR, PC_guess_LR, PI_LR, w_LR, laboralloc_LR, p_KR_bar_LR, KR_LR, KF_LR, p_KR_LR_S, p_KR_LR_W, p_B,
    p_F_LR, Lsector, YF_LR, curtailmentfactor, diffK, diffp, KRshifter, P_out, result_Dout_LR, result_Yout_LR, result_YFout_LR,
    result_price_LR, KP_LR = solve_power_output(RWParams, params, RunBatteries, RunCurtailment,
                                                Initialprod, R_LR, majorregions, Linecounts, linconscount,
                                                regionParams, curtailmentswitch)



# Get fossil fuel usage in the initial steady state
#YF_LR = max.(YE_LR .- regionParams.theta .* KR_LR, 0)
FUtilization = YF_LR ./ KF_LR
#mean(FUtilization[1:majorregions.n[1]]))

# Compute electricity and fossil fuel usage in industry and electricity sectors
e2_LR = laboralloc .* repeat(D_LR, 1, params.I) .* (repeat(params.Vs[:,2]', params.J, 1) ./ (repeat(params.Vs[:,2]', params.J, 1) + repeat(params.Vs[:,3]', params.J, 1)))
fusage_ind_LR = sum(e2_LR, dims=2) .* (p_E_LR / p_F_LR) .^ params.psi
fusage_power_LR = (YF_LR ./ KF_LR .^ params.alpha2) .^ (1 / params.alpha1)
fusage_total_LR = sum(fusage_power_LR) + sum(fusage_ind_LR)

# compute fossil usage as a share of GDP
GDP = sum(w_LR .* params.L .+ p_E_LR .* D_LR .+ rP_LR .* KP_LR .+ p_F_LR .* fusage_ind_LR)
Fossshare = sum(p_F_LR .* fusage_ind_LR) / GDP

# Save long run variables for transition
I_tota_S_LR = copy(Depreciation_LR_S)
I_tota_W_LR = copy(Depreciation_LR_W)

if hoursofstorage == 0
    @save "Guesses/laboralloc_LR_guess.jld2" laboralloc_LR
    @save "Guesses/KR_LR_S_guess.jld2" KR_LR_S
    @save "Guesses/KR_LR_W_guess.jld2" KR_LR_W
    @save "Guesses/p_E_LR_guess.jld2" p_E_LR
    @save "Guesses/w_LR_guess.jld2" w_LR
    @save "Guesses/Dout_guess_LR.jld2" result_Dout_LR
    @save "Guesses/Yout_guess_LR.jld2" result_Yout_LR
    @save "Guesses/PC_guess_LR.jld2" PC_guess_LR
end


wageresults[:,2] = w_real
wr = (wageresults[:, 2] ./ wageresults[:, 1]) .- 1
wagechange = [wr regions.csr_id]

priceresults_LR = p_E_LR

# Get changes in welfare from the different components
welfare_wagechange = (log.(w_LR ./ PC_guess_LR) .- log.(wage_init ./ PC_guess_init)) .* (wage_init .* params.L ./ Expenditure_init)
welfare_capitalchange = (log.(KP_LR ./ PC_guess_LR) .- log.(KP_init ./ PC_guess_init)) .* (rP_init .* KP_init ./ Expenditure_init)

welfare_electricitychange = (log.((R_LR .* KR_LR .* p_KR_bar_LR .* PC_LR + R_LR .* KF_LR .* PC_LR) ./ PC_LR) .- 
    log.((R_LR .* (KR_init_W + KR_init_S) .* p_KR_bar_init .* PC_init + R_LR .* KF_init .* PC_init))) .* 
    ((1 - params.beta) * (R_LR .* (KR_init_W + KR_init_S) .* p_KR_bar_init .* PC_init + R_LR .* KF_init .* PC_init) ./ Expenditure_init)

welfare_fossilchange = -fossilsales ./ Expenditure_init

# Get changes in manufacturing share
comparativeadvantagechange = laboralloc_init ./ laboralloc_LR
scatter(log.(params.L[majorregions.rowid2[1]:majorregions.rowid[1]]), result_price_LR[1])

w_real_LR = w_LR ./ PC_guess_LR
GDP_ind_LR = (w_LR .* params.L .+ p_E_LR .* D_LR .+ rP_LR .* KP_LR .+ p_F_LR .* fusage_ind_LR) ./ PC_guess_LR

end