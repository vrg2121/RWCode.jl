# Market.jl
module Market
using Main.DataAdjustments, Main.Params, MAT, Main.ParameterizeRun, Main.DataLoads, Printf, Random, LinearAlgebra, SparseArrays, JuMP, Ipopt, Distributed, JLD2, Statistics, DataFrames

import Main.DataLoads: wage_init, regionParams, KR_init_S, KR_init_W, R_LR, sectoralempshares
import Main.Params: params, majorregions, Linecounts, kappa, regions, linconscount, updw_w, upw_z

import Main.MarketEquilibrium: mycon, obj, Price_Solve, wage_update_ms, obj2, calc_subsidyUS
using Main.MarketFunctions

export result_Yout_init, laboralloc, wageresults, PC_guess_init, Expenditure_init, 
	   KP_init, rP_init, p_KR_bar_init, PC_init, KF_init, fossilsales, laboralloc_init, 
	   p_E_init, D_init, YE_init, PI_init, p_KR_init_S, p_KR_init_W, p_F_path_guess, 
	   wedge, priceshifterupdate, p_F, renewshareUS, p_F_int


function solve_initial_equilibrium(params, wage_init::Vector{Float64}, majorregions::DataFrame,
        regionParams, KR_init_S::Matrix{Float64}, KR_init_W::Matrix{Float64}, R_LR, sectoralempshares::Matrix,
        Linecounts::DataFrame, kappa::Int, regions::DataFrame, linconscount::Int, updw_w, upw_z)
        
	# allocations for variables
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

	# load variable guesses
	w_guess, p_E_init, result_Pout_init, laboralloc, 
	PC_guess_init, result_Dout_init, result_Yout_init, 
	p_F_path_guess, wedge, priceshifterupdate, fossilsales = open_mat_var()


	# initial fossil sales guess
	p_F = 0.05

	Lsector = laboralloc.*params.L

	diffend=1

        
    while diffend>10^(-2)
		YF_init = Vector{Float64}(undef, 2531)
		PI_init = Vector{Float64}(undef, 2531)
		KP_init = Vector{Float64}(undef, 2531)
		

		# optimization for each reagion kk
		rP_init, pg_init_s, pE_market_init = market_setup!(rP_init, pg_init_s, pE_market_init, wage_init, params, p_E_init, p_F, R_LR, PC_guess_init)

        #Threads.@threads for kk in 1:(params.N - 1)
            kk = 1
			# the local call in this section localizes the memory to each thread to reduce crossing data 
			
			# set up optimization problem for region kk
			ind = majorregions.rowid2[kk]:majorregions.rowid[kk]
			n = majorregions.n[kk]
			l_ind = Linecounts.rowid2[kk]:Linecounts.rowid[kk]
			gam = RWParams.Gam[kk]
			l_n = Linecounts.n[kk]

			@views secalloc = laboralloc[ind, :]
			@views Lshifter = Lsector[ind, :]
			Kshifter = Lsector[ind, :] .* (params.Vs[:,4]' .* ones(n, 1)) ./
					(params.Vs[:,1]' .* ones(n, 1)) .*
					(wage_init[ind] ./
					rP_init[ind])
			#Ltotal=sum(Lshifter,dims=2)
			#Jlength= n
			# all of these are varying correctly

			# define data for inequality constraints
			linecons = copy(RWParams.Zmax[l_ind])
			Gammatrix = hcat(zeros(size(gam, 1)), gam)
			
			if linconscount < l_n
					Random.seed!(1)  
					randvec = rand(l_n)  
					randvec = randvec .> (linconscount / l_n)  
					Gammatrix[findall(randvec), :] .= 0 
			end   
			
			id_mx = Matrix{Float64}(I, n, n)
			stacker = hcat(-id_mx, id_mx)
			Gammatrix = sparse(Gammatrix * stacker)
			Gammatrix = vcat(Gammatrix, -Gammatrix)
			linecons = vcat(linecons, linecons)
			# all of these are calculated correctly for each kk

			# define shifters for objective function
			@views pg_s = pg_init_s[ind, :]
			@views prices = (pE_market_init[ind] .* ones(1, params.I))
			@views power = (params.Vs[:,2]' .* ones(n, 1)) + (params.Vs[:,3]' .* ones(n, 1))


			@views p_1 = pg_s .* (kappa .+ (prices ./ (kappa .* p_F)).^(params.psi - 1)).^(params.psi ./ (params.psi - 1) .* (params.Vs[:, 3]' .* ones(n, 1)))
			@views p_2 = (1 .+ (params.Vs[:,3]' .* ones(n, 1)) ./ (params.Vs[:, 2]' .* ones(n, 1))).^(-(params.Vs[:,2]' .* ones(n, 1)) .- (params.Vs[:,2]' .* ones(n, 1)))
			@views p_3 = params.Z[ind] .* params.zsector[ind, :]
			@views p_4 = Lshifter.^(params.Vs[:,1]' .* ones(n, 1)) .* Kshifter .^ (params.Vs[:,4]' .* ones(n, 1))
			shifter = @. p_1 * p_2 * p_3 * p_4
			shifter = @. shifter * secalloc ^ power 

			@views KRshifter = @. regionParams.thetaS[ind] * KR_init_S[ind] + 
					regionParams.thetaW[ind] * KR_init_W[ind]
			@views KFshifter=regionParams.KF[ind]

			# define bounds
			YFmax=regionParams.maxF[ind]
			LB = [zeros(n); KRshifter]
			UB = [fill(1000, n); YFmax + KRshifter]

			guess = [KRshifter .- 0.001; KRshifter]
			l_guess = length(guess)

			# solve market equilibrium
			## constraint and objective functions are in MarketEquilibrium.jl
			x = Vector{Float64}(undef, l_guess)
            mid = l_guess ÷ 2
			model = Model(Ipopt.Optimizer)
			set_silent(model)
			@variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
           
            @expression(model, Dvec, x[1:mid])
            @expression(model, Yvec, x[mid+1:end])
            @expression(model, Pvec, Yvec .- Dvec)
            @expression(model, sPvec, sum(Pvec))
            @expression(model, quad_mat, -params.Rweight .* Pvec[2:end]' * regionParams.B[kk] * Pvec[2:end])
            add_to_expression!(quad_mat[1], sPvec-Pvec[1]^2)

            @constraint(model, c, quad_mat[1] == 0)
            
			#@constraint(model, c1, mycon(x, regionParams.B[kk], params) .== 0)
			@objective(model, Min, obj(x, power, shifter, KFshifter, KRshifter, p_F, params))
			optimize!(model)

			## save output as P_out
			P_out = value.(x)

			result_price_init[kk].=Price_Solve(P_out, shifter, n, params) #.MarketEquilibrium.jl
			mid = length(P_out) ÷ 2
			@views result_Dout_init[kk] .= P_out[1:mid]
			@views result_Yout_init[kk] .= P_out[1+mid:end]
			Pvec = P_out[1+mid:end] .- P_out[1:mid]
			Losses = Pvec[2:end]' * regionParams.B[kk] * Pvec[2:end] .* params.Rweight
			@views Lossfac_init[kk] = Losses ./ sum(result_Yout_init[kk])


        #end   


		# solve market equilibrium
		kk = params.N

		KRshifter, YFmax, guess, power, KFshifter, kk, shifter = second_loop(kk, majorregions, laboralloc, Lsector, params, wage_init, rP_init,
																				result_Pout_init, pg_init_s, pE_market_init, kappa, 
																				regionParams, KR_init_S, KR_init_W, p_F)


        Threads.@threads for jj in 1:majorregions.n[kk]
                local con = [1 -1]
                local guess = [1; KRshifter[jj]]
                local LB = [0; KRshifter[jj]]
                local UB = [10^6; YFmax[jj] + KRshifter[jj]]
                local l_guess = length(guess)

				local x = Vector{Float64}(undef, 2)
                model = Model(Ipopt.Optimizer)
                set_silent(model)
                @variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
                @constraint(model, c1, con * x <= 0) 
                @objective(model, Min, obj2(x, power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F, params))
                optimize!(model)

				local P_out1 = value.(x)  # Default to local scope

				if jj == majorregions.n[kk]
					global P_out = P_out1  # Only set to global on specific condition
				end 

                local result_Dout_init[kk][jj] = P_out1[1]
                local result_Yout_init[kk][jj] = P_out1[2]

                result_price_init[kk][1][1, jj] = Price_Solve(P_out1, shifter[jj], 1, params)[1]
				
				# put bounds on prices
                result_price_init[kk][1][1, jj] = max(result_price_init[kk][1][1, jj], 0.001)
                result_price_init[kk][1][1, jj] = min(result_price_init[kk][1][1, jj], 0.3)
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
					matrix_values = result_price_init[kk][1]
					p_E_init[ind] .= vec(matrix_values)
					@views D_init[ind].= vec(result_Dout_init[kk])
					@views YE_init[ind] .= vec(result_Yout_init[kk])
					@views PI_init[ind, 1] .= (sum(matrix_values.*(result_Dout_init[kk] .- result_Yout_init[kk]))) .*
										(params.L[ind, 1]) ./
										(sum(params.L[ind, 1]))
			end

        end


        #global PI_init = map(y -> let x = round(y; digits=6); x == -0.0 ? 0.0 : x end, PI_init)
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
        #diffsec = maximum(abs.(sectoralempshares[:, 2:end] .- laboralloc[:, 1:end]), dims=1)
        #diff_w2 = maximum(abs.(w_update .- w_guess) ./ w_update)

        # update relative labor allocations
        laboralloc .= Lsector ./ params.L
        relab .= laboralloc ./ sum(laboralloc, dims=2)
		relexp .= Xjdashs ./ sum(Xjdashs, dims = 2)
        Lsector .= Lsector .* clamp.(1 .+ 0.05 .* (relexp .- relab), 0.9, 1.1)
        Lsector .= Lsector ./ sum(Lsector, dims=2) .* params.L

        # calibrate sectoral params.Z to move sectoral emp closer
        params.zsector[:, :] .= params.zsector[:, :] .* clamp.(1 .+ upw_z .* (sectoralempshares[:, 2:end] .- laboralloc[:, :]), 0.99, 1.01)

        # Update consumption price guess
        PC_guess_init .= 0.2 .* PC_init .+ (1 - 0.2) .* PC_guess_init

        # get fossil fuel usage in the initial steady state
		global @. YF_init = YE_init - regionParams.thetaS * KR_init_S - regionParams.thetaW * KR_init_W
        FUtilization .= YF_init ./ regionParams.maxF
        #mean(FUtilization[1:majorregions[1, :n]])

        # calibrate efficiency wedges costs to match initial power prices across countries
		# set subsidy amount for transition
		subsidy_US = calc_subsidyUS(p_E_init, regions, majorregions) 


		fossilsales, Expenditure_init = elec_fuel_expenditure!(laboralloc, D_init, params, p_E_init, p_F, YF_init, regionParams,
																		wage_init, rP_init, KP_init, fossilsales, YE_init, regions, PI_init)
		println("Initial Calibration= ", diffend)

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
			p_F = p_F)
end

mrkteq = solve_initial_equilibrium(params, wage_init, majorregions, regionParams, KR_init_S,
                                KR_init_W, R_LR, sectoralempshares, Linecounts, kappa, regions, 
                                linconscount, updw_w, upw_z);


Z = params.Z
Zsec = params.zsector
w_guess = mrkteq.w_guess
p_E_init = mrkteq.p_E_init
result_Dout_init = mrkteq.result_Dout_init
result_Yout_init = mrkteq.result_Yout_init
PC_guess_init = mrkteq.PC_guess_init
laboralloc = mrkteq.laboralloc
wedge = mrkteq.wedge
priceshifterupdate = mrkteq.priceshifterupdate
fossilsales = mrkteq.fossilsales


# Save variables
@save "Guesses/w_guess_mat.jld2" w_guess
@save "Guesses/p_E_guessmat.jld2" p_E_init
@save "Guesses/Dout_guess_init.jld2" result_Dout_init
@save "Guesses/Yout_guess_init.jld2" result_Yout_init
@save "Guesses/PC_guess_init.jld2" PC_guess_init

@save "Guesses/laboralloc_guess.jld2" laboralloc
@save "Guesses/z_mat.jld2" Z
@save "Guesses/z_sec_mat.jld2" Zsec
@save "Guesses/wedge_vec.jld2" wedge
@save "Guesses/priceshifterupdate_vec.jld2" priceshifterupdate
@save "Guesses/fossilsales_guess.jld2" fossilsales

# set initial power output vector
P_out_init = mrkteq.P_out

# set Q initial to be the solar already installed
Qtotal_init_S = sum(KR_init_S)
Qtotal_init_W = sum(KR_init_W)
p_KR_init_S=(Initialprod+Qtotal_init_S).^(-params.gammaS);    
p_KR_init_W=(Initialprod+Qtotal_init_W).^(-params.gammaW);
SShare_init=(regionParams.thetaS./p_KR_init_S).^params.varrho./((regionParams.thetaS./p_KR_init_S).^params.varrho+(regionParams.thetaW./p_KR_init_W).^params.varrho)
thetabar_init = regionParams.thetaS .* SShare_init + regionParams.thetaW .* (1 .- SShare_init)
p_KR_bar_init = SShare_init .* p_KR_init_S + (1 .- SShare_init) .* p_KR_init_W
pE_FE_init=(p_KR_bar_init-p_KR_bar_init*(1-params.deltaR)./R_LR).*regionParams.costshifter./thetabar_init

wageresults = Matrix{Float64}(undef, 2531, 2)
priceresults = Vector{Float64}(undef, 2531)
PCresults = Vector{Float64}(undef, 2531)
KF_init = Vector{Float64}(undef, 2531)

wageresults[:,1].=copy(mrkteq.W_Real)
priceresults[:,1].=copy(p_E_init)
laboralloc_init=copy(laboralloc)
PCresults[:,1].=copy(p_E_init)


#get renewable shares
renewshareUS=1-(sum(YF_init[1:majorregions.n[1],:])./sum(mrkteq.YE_init[1:majorregions.n[1]]));
renewshareEU=1-(sum(YF_init[majorregions.rowid[1]+1:majorregions.rowid[2]]))./sum(mrkteq.YE_init[majorregions.rowid[1]+1:majorregions.rowid[2]])

p_F_int=copy(mrkteq.p_F)
KF_init=copy(regionParams.KF)
end