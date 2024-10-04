module MarketEquilibrium
export hessinterior, obj, mycon, obj2, Price_Solve, wage_update_ms, grad_f

import ..ParamsFunctions: StructParams


function obj(Inputvec::Vector, power::Matrix, shifter::Matrix, KFshifter::Union{SubArray, Vector}, 
        KRshifter::Vector, p_F::Union{Float64, Vector, Int}, params::StructParams)
    mid = length(Inputvec) ÷ 2
    Dvec = Inputvec[1:mid]
    Yvec = Inputvec[1+mid:end]
    #Dsec = Dvec .* ones(1, params.I)
    Dsec = broadcast(*, Dvec, ones(1, params.I))

    power2 = 1 / params.alpha1
    value1 = -(Dsec .^ power) .* shifter
    value1 = sum(value1, dims=2)
    value2 = @. p_F * ((Yvec - KRshifter) / KFshifter ^ params.alpha2) ^ power2
    value = sum(value1) + sum(value2)    
    return value
    
end

function obj2(Inputvec::Vector, power::Float64, shifter::Float64, KFshifter::Float64, 
            KRshifter::Float64, p_F::Union{Float64, Int64}, params::StructParams)
    mid = length(Inputvec) ÷ 2
    Dvec = Inputvec[1:mid]
    Yvec = Inputvec[1+mid:end]
    Dsec = repeat(Dvec, 1, params.I)
    power2 = 1 / params.alpha1
    value1 = -Dsec .^ power .* shifter
    value2 = @. p_F * ((Yvec - KRshifter) / KFshifter ^ params.alpha2) ^ power2
    value1 = sum(value1, dims=2)
    value = sum(value1) + sum(value2)

    return value
    
end

function Price_Solve(Inputvec::Vector{Float64}, shifter::Union{Matrix, Float64}, Jlength::Int64, params::StructParams)
    mid = length(Inputvec) ÷ 2
    Dvec = Inputvec[1:mid]

    prices = sum(((params.Vs[:,2]' .* ones(Jlength, 1)) .+ 
            (params.Vs[:,3]' .* ones(Jlength, 1))) .* ((Dvec) .^ ((params.Vs[:,2]' .* ones(Jlength, 1)) .+ 
            (params.Vs[:,3]' .* ones(Jlength, 1)) .- 1)) .* shifter, dims=2)

    return vec(prices)
end


function location_prices!(pijs::Vector{Matrix{Float64}}, PCs::Matrix{Float64}, Xjdashs::Matrix{Float64}, Yjdashs::Matrix{Float64},
    w0::Vector{Float64}, p_E_D::Vector{Float64}, params::StructParams, Ej::Matrix{Float64}, 
    p_F::Union{Float64, Vector, Int}, r::Matrix{Float64})

    for i = 1:params.I
        ttau = params.tau[i]

        @views tpijs = ttau .* w0 .^ params.Vs[i, 1] .* p_E_D .^ (params.Vs[i, 2] + params.Vs[i, 3]) .* 
        (params.kappa + (params.kappa .* p_F ./ p_E_D) .^ (1 - params.psi)) .^ (-(params.psi ./ (params.psi - 1)) .* params.Vs[i, 3]) .*
        r .^ params.Vs[i, 4] ./
        (params.Z .* params.zsector[:, i] .* params.cdc)

        PCs[:, i] = (sum(tpijs .^ (1 - params.sig), dims=1)) .^ (1 / (1 - params.sig))
        @views Xjdashs[:, i] = sum((tpijs) .^ (1 - params.sig) .* (params.betaS[:, i] .* Ej ./ PCs[:, i].^(1 - params.sig))', dims=2)
        @views Yjdashs[:, i] = sum((tpijs).^(-params.sig) .* (params.betaS[:, i] .* Ej ./ PCs[:, i].^(1 - params.sig))', dims=2)
        pijs[i] .= tpijs
    end
end


function update_wage_data!(tpijs::Matrix, params::StructParams, w0::Union{Matrix, Vector}, pES::Union{Vector, Matrix},
    p_F::Union{Float64, Vector, Int}, r::Union{Vector, Matrix}, Ej, PCs::Matrix, Xjdashs::Matrix,
    Yjdashs::Matrix)

    for i = 1:params.I
        @views tpijs .= params.tau[i] .* 
                        w0 .^ params.Vs[i, 1] .* 
                        pES .^ (params.Vs[i, 2] + params.Vs[i, 3]) .* 
                        (params.kappa + (params.kappa .* p_F ./ pES) .^ (1 - params.psi)) .^ 
                        (-(params.psi / (params.psi - 1)) * params.Vs[i, 3]) .*
                        r .^ params.Vs[i, 4] ./
                        (params.Z .* params.zsector[:, i] .* params.cdc)

        PCs[:, i] = (sum(tpijs .^ (1 - params.sig), dims=1)) .^ (1 / (1 - params.sig))

        factor = (params.betaS[:, i] .* Ej ./ PCs[:, i].^(1 - params.sig))'
        Xjdashs[:, i] = sum(tpijs .^ (1 - params.sig) .* factor, dims=2)
        Yjdashs[:, i] = sum(tpijs .^ (-params.sig) .* factor, dims=2)
    end

end

function price_adjustments!(PC::Vector{Float64}, PCs::Array{Float64}, params::StructParams, w0::Union{Vector, Matrix}, 
                        Xjdashs::Array{Float64}, Xj::Union{Vector, Matrix}, pES::Union{Vector, Matrix}, 
                        pED::Union{Vector, Matrix}, p_F::Union{Float64, Vector, Int}, W_Real::Vector{Float64}, 
                        w_adjustment_factor::Union{Vector, Matrix}, Xjdash::Matrix{Float64})
    updw=0.5
    PC .= prod(PCs .^ params.betaS, dims=2)
    Xjdash .= sum(Xjdashs, dims=2)
    w_adjustment_factor .= w0 .* min.(max.(1 .+ updw .* (Xjdash .- Xj) ./ Xj, 0.2), 1.1)
    w0 .= w_adjustment_factor ./ (w_adjustment_factor[1])
    pES .= pES ./ (w_adjustment_factor[1])
    pED .= pED ./ (w_adjustment_factor[1])
    p_F = p_F / (w_adjustment_factor[1])
    W_Real .= w0 ./ PC
end

function wage_update_ms(w::Union{Vector, Matrix}, p_E_D::Union{Vector, Matrix},p_E_S::Union{Vector, Matrix}, p_F::Union{Float64, Vector, Int}, D_E::Vector, Y_E::Vector, 
    r::Union{Vector, Matrix}, KP::Union{Vector, Matrix}, Pi::Vector, fossil::Union{Vector, Matrix, Int}, params::StructParams)
    
    PCs = Array{Float64}(undef, 2531, 10)
    Xjdashs = Array{Float64}(undef, 2531, 10)
    Yjdashs = Array{Float64}(undef, 2531, 10)
    w_adjustment_factor = Matrix{Float64}(undef, 2531, 1)
    W_Real = Vector{Float64}(undef, 2531)
    PC = Vector{Float64}(undef, 2531)
    tpijs = Matrix{Float64}(undef, 2531, 2531)
    Xjdash = Matrix{Float64}(undef, 2531, 1)

    w0 = copy(w)
    pED = copy(p_E_D)
    pES = copy(p_E_S)

    Xj = @. w0 * params.L + pED * D_E + r * KP  
    Ej = @. w0 * params.L + pES * Y_E + r * KP + Pi + fossil

    update_wage_data!(tpijs, params, w0, pED, p_F, r, Ej, PCs, Xjdashs, Yjdashs)
    
    # calculate price indices and adjust wages
    price_adjustments!(PC, PCs, params, w0, Xjdashs, Xj, pES, pED, p_F, W_Real, w_adjustment_factor, Xjdash)

    return w0, W_Real, sum(Xj), PC, Xjdashs, PCs, Yjdashs, Xj

end


end