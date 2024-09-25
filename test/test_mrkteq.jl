module test_mrkteq

function mycon(Inputvec::Vector, mat::Matrix, params)
    mid = length(Inputvec) ÷ 2

    Dvec = @view Inputvec[1:mid]
    Yvec = @view Inputvec[1+mid:end]
    Pvec = Yvec .- Dvec

    #@views ceq = sum(Pvec) - params.Rweight * (Pvec[2:end]' * mat * Pvec[2:end]) - Pvec[1]^2
    return sum(Pvec) - params.Rweight * (Pvec[2:end]' * mat * Pvec[2:end]) - Pvec[1]^2
end



function obj(Inputvec::Vector, power::Matrix, shifter::Matrix, KFshifter::Union{SubArray, Vector}, 
        KRshifter::Vector, p_F, params)
    mid = length(Inputvec) ÷ 2
    Dvec = Inputvec[1:mid]
    Yvec = Inputvec[1+mid:end]
    Dsec = Dvec .* ones(1, params.I)

    power2 = 1 / params.alpha1
    value1 = -(Dsec .^ power) .* shifter
    value2 = @. p_F * ((Yvec - KRshifter) / KFshifter ^ params.alpha2) ^ power2
    value1 = sum(value1, dims=2)
    #value = sum(value1) + sum(value2)

    return sum(value1) + sum(value2)
    
end

function main_model()
    local x = Vector{Float64}(undef, l_guess)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_string_names_on_creation(model, false)

    @variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
    
    @constraint(model, c1, mycon(x, regionParams.B[kk], params) == 0) # 4156 calls on profiler, 1504 calls to mycon, marketequilibrium 19
    @objective(model, Min, obj(x, power, shifter, KFshifter, KRshifter, p_F_in, params)) # 374 calls on profiler
    optimize!(model)

    return value.(x)
end

end