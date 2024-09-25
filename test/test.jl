s = 21
larry = 13

function larry_1()
    global larry1 = 13
    while larry1 < 20
        larry1+=1
    end 
    return larry1
    # global variables must be noted as global within the function
    # return a variable does not bring it into global scope unless specified as global
end

larry1 = larry_1()

function larry_test()
while s > 15
    b = 5
    for kk = 5:9
        global a = Vector{Float64}(undef, kk)
        a .= 15-kk
    end
    a .= b - 1 + larry
    a .-= 2
    s = sum(a)
    global larry -= 1
end
end
larry_test()
println(a)


lsec = collect(range(50, stop=60, length=10))


kk = 1
for kk = 1
    lsec .-=2
end



kk = 1
for kk = 1
    tsec = collect(range(50, stop=60, length=10))
    tsec .-=2
end


z_file = matopen("Guesses/z_mat.mat")
z_check = read(z_file, "Z")
close(z_file)



t = 100
rosenbrock(t, x...) = (1 - x[1])^2 + t * (x[2] - x[1]^2)^2
function ∇rosenbrock(g::AbstractVector, x...)
    g[1] = 400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2
    g[2] = 200 * (x[2] - x[1]^2)
    return
end

# testing threading in market
# MULTI THREADING WORKS HERE
pg_init_s .= wage_init .^ repeat(params.Vs[:, 1]', outer=[params.J, 1]) .* 
                p_E_init .^ (repeat(params.Vs[:, 2]', outer=[params.J, 1]) + repeat(params.Vs[:, 3]', outer=[params.J, 1])) .* 
                (params.kappa + (params.kappa .* p_F ./ p_E_init) .^ (1-params.psi)) .^ (-(params.psi / (params.psi - 1)) * params.Vs[:, 3]') .* 
                rP_init .^ repeat(params.Vs[:, 4]', outer=[params.J, 1]) ./ 
                (params.Z .* params.zsector .* params.cdc)
        
pE_market_init .= copy(p_E_init)


using Base.Threads

lk = ReentrantLock()
jj=0
tol = 1
while tol < 3
    global lk
    while jj < 2
        Threads.@threads for kk in 1:(params.N - 1) 
            #for kk = 1:(params.N - 1)
            # to speed up calculation, pre-allocate data slices in loop
            #nr = majorregions.rowid[kk] - majorregions.rowid2[kk] + 1
            #Kshifter = Matrix{Float64}(undef, nr, size(Lsector, 2))

            # set up optimization problem for region kk
            secalloc = laboralloc[majorregions.rowid2[kk]:majorregions.rowid[kk], :]
            Lshifter = Lsector[majorregions.rowid2[kk]:majorregions.rowid[kk], :]
            Kshifter = Lsector[majorregions.rowid2[kk]:majorregions.rowid[kk], :] .* repeat(params.Vs[:,4]', majorregions.n[kk], 1) ./
                    (repeat(params.Vs[:,1]', majorregions.n[kk], 1)) .*
                    (wage_init[majorregions.rowid2[kk]:majorregions.rowid[kk]] ./
                    rP_init[majorregions.rowid2[kk]:majorregions.rowid[kk]])
            Ltotal=sum(Lshifter,dims=2)
            Jlength= majorregions.n[kk]


            # define data for inequality constraints
            #nr1 = Linecounts.rowid[kk] - Linecounts.rowid2[kk] + 1
            #linecons = Vector{Float64}(undef, nr1)

            linecons = copy(RWParams.Zmax[Linecounts.rowid2[kk]:Linecounts.rowid[kk]])
            Gammatrix = hcat(zeros(size(RWParams.Gam[kk], 1)), RWParams.Gam[kk])

            if linconscount < Linecounts.n[kk]
                    Random.seed!(1)  
                    randvec = rand(Linecounts.n[kk])  
                    randvec = randvec .> (linconscount / Linecounts.n[kk])  
                    Gammatrix[findall(randvec), :] .= 0 
            end   

            n = majorregions.n[kk]
            id_mx = Matrix{Float64}(I, n, n)
            stacker = hcat(-id_mx, id_mx)

            Gammatrix = sparse(Gammatrix * stacker)
            Gammatrix = vcat(Gammatrix, -Gammatrix)
            linecons = vcat(linecons, linecons)

            # define shifters for objective function
            pg_s = pg_init_s[majorregions.rowid2[kk]:majorregions.rowid[kk], :]
            prices = repeat(pE_market_init[majorregions.rowid2[kk]:majorregions.rowid[kk]], 1, params.I)
            power = repeat(params.Vs[:,2]',Jlength, 1 ) + repeat(params.Vs[:,3]', Jlength,1)

            p_1 = pg_s .* (kappa .+ (prices ./ (kappa .* p_F)).^(params.psi .- 1)).^(params.psi ./ (params.psi .- 1) .* repeat(params.Vs[:, 3]', majorregions.n[kk], 1))
            p_2 = (1 .+ repeat(params.Vs[:, 3]', outer = (majorregions.n[kk], 1)) ./ repeat(params.Vs[:, 2]', outer = (majorregions.n[kk], 1))).^(-repeat(params.Vs[:, 2]', outer = (majorregions.n[kk], 1)) .- repeat(params.Vs[:, 2]', outer = (majorregions.n[kk], 1)))
            p_3 = params.Z[majorregions.rowid2[kk]:majorregions.rowid[kk]] .* params.zsector[majorregions.rowid2[kk]:majorregions.rowid[kk], :]
            p_4 = Lshifter.^repeat(params.Vs[:, 1]', outer = (majorregions.n[kk], 1)) .* Kshifter.^repeat(params.Vs[:, 4]', outer = (majorregions.n[kk], 1))
            shifter = p_1 .* p_2 .* p_3 .* p_4
            shifter = shifter .* secalloc .^ power 

            KRshifter = regionParams.thetaS[majorregions.rowid2[kk]:majorregions.rowid[kk]] .* KR_init_S[majorregions.rowid2[kk]:majorregions.rowid[kk]] .+ 
                    regionParams.thetaW[majorregions.rowid2[kk]:majorregions.rowid[kk]] .* KR_init_W[majorregions.rowid2[kk]:majorregions.rowid[kk]]
            KFshifter=regionParams.KF[majorregions.rowid2[kk]:majorregions.rowid[kk]]

            # define bounds
            YFmax=regionParams.maxF[majorregions.rowid2[kk]:majorregions.rowid[kk]]
            LB = [zeros(majorregions.n[kk]); KRshifter]
            UB = [fill(1000, majorregions.n[kk]); YFmax + KRshifter]

            guess = [KRshifter .- 0.001; KRshifter]
            l_guess = length(guess)

            # solve market equilibrium
            ## constraint and objective functions are in MarketEquilibrium.jl
            x = []
            model = Model(Ipopt.Optimizer)
            set_silent(model)
            @variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
            @constraint(model, c1, mycon(x, regionParams.B[kk], params) == 0)
            @objective(model, Min, obj(x, power, shifter, KFshifter, KRshifter, p_F, params))
            optimize!(model)
            ## save output as P_out
            P_out = value.(x)

            lock(lk) do
                result_price_init[kk]=Price_Solve(P_out, shifter, Jlength, params) #.MarketEquilibrium.jl
            end

            mid = length(P_out) ÷ 2

            lock(lk) do
                result_Dout_init[kk] .= P_out[1:mid]
            end

            lock(lk) do
                result_Yout_init[kk] .= P_out[1+mid:end]
            end

            Pvec = P_out[1+mid:end] .- P_out[1:mid]
            Losses = Pvec[2:end]' * regionParams.B[kk] * Pvec[2:end] .* params.Rweight

            lock(lk) do
                Lossfac_init[kk] = Losses / sum(result_Yout_init[kk])
            end
        end
        jj += 1     # loop should run twice
        println("Inner while loop complete: jj = $jj")
    end
    tol += 1
    println("Outer while loop complete: tol = $tol")
end

# threading checked
using Base.Threads

function parallel_processing_with_condition(max_value)
    continue_processing = true
    current_max = 0

    # Initialize an array to store results of parallel computation with an initial size
    results = Array{Int}(undef, 10)  # Assuming the batch size is always 10

    while continue_processing
        # Generate a batch of random numbers for this iteration
        data_batch = rand(1:100, 10)

        # Use threading to process each item in the batch in parallel
        Threads.@threads for i in 1:length(data_batch)
            # Simulate some processing by squaring the number
            results[i] = data_batch[i]^2
        end

        # Check results to decide if we should continue
        current_max = maximum(results)
        println("Current batch max squared value: $current_max")
        continue_processing = current_max < max_value

        # Here, we assume some condition based on results or other factors
        if !continue_processing
            println("Maximum value reached or exceeded: $current_max")
        else
            println("Continuing processing...")
        end
    end
end

# Example usage:
parallel_processing_with_condition(5000)




# dummy version of market.jl multithreading

d1 = Matrix{Float64}(undef, 10, 2)
d2 = Matrix{Float64}(undef, 10, 2)
scalar = 3

t = 100
r = 400
s = 200
f(t, x...) = (1 - x[1])^2 + t * (x[2] - x[1]^2)^2
function ∇f(g, r, s, x...)
    g[1] = r * x[1]^3 - r * x[1] * x[2] + 2 * x[1] - 2
    g[2] = s * (x[2] - x[1]^2)
    return
end
function ∇²f(H, r, s, x...)
    H[1, 1] = 1200 * x[1]^2 - r * x[2] + 2
    # H[1, 2] = -400 * x[1]  <-- Not needed. Fill the lower-triangular only.
    H[2, 1] = -r * x[1]
    H[2, 2] = s
    return
end

model = Model(Ipopt.Optimizer)
register(model, :rosenbrock, 2, f, ∇f, ∇²f)
@variable(model, x[1:2])
@NLobjective(model, Min, rosenbrock(t, x[1], x[2]))
optimize!(model)


######
# simplified working model
using Random
using SparseArrays
function obj_f(power, shifter, KFshifter, KRshifter, p_F, params, Inputvec)
    mid = length(Inputvec) ÷ 2
    Dvec = Inputvec[1:mid]
    Yvec = Inputvec[1+mid:end]
    Dsec = repeat(Dvec, 1, params.I)
    power2 = 1 / params.alpha1
    value1 = -Dsec .^ power .* shifter
    value2 = p_F .* ((Yvec .- KRshifter) ./ KFshifter .^ params.alpha2) .^ power2
    value1 = sum(value1, dims=2)
    value = sum(value1) + sum(value2)

    return value
end

function grad_f(power, shifter, KFshifter, KRshifter, p_F, params, Inputvec)
    mid = length(Inputvec) ÷ 2
    Dvec = Inputvec[1:mid]
    Yvec = Inputvec[1+mid:end]
    Dsec = repeat(Dvec, 1, params.I)
    power2 = 1 / params.alpha1

    grad1f = -sum((power).* Dsec .^ (power .- 1) .* shifter, dims=2)
    grad2f = power2 .* p_F .* (Yvec .- KRshifter) .^ (power2-1) .* (1 ./ KFshifter .^ params.alpha2).^ (power2)
    gradf = sparse([grad1f, grad2f])
    return gradf
end

function f_hess(lambda, power, shifter, KFshifter, KRshifter, p_F, mat, params, Inputvec)
    Dvec = Inputvec[1:end÷2]  # ÷ is integer division
    Yvec = Inputvec[end÷2+1:end]
    J = length(Dvec)

    Dsec = repeat(Dvec, 1, params.I)
    power2 = (1 / params.alpha1)

    piece1 = -sum(((power .- 1) .* power) .* Dsec .^ (power .- 2) .* shifter, dims=2)
    piece1 = Diagonal(piece1)

    piece2 = (power2 - 1) .* power2 .* p_F .* (Yvec .- KRshifter) .^ (power2 - 2) .* (1 ./ KFshifter .^ params.alpha2) .^ power2
    piece2 = Diagonal(piece2)

    f2 = [piece1 zeros(J, J); zeros(J, J) piece2]
    f2 = sparse(f2)

    xmat = params.Rweight * 2 * mat
    zeros1block = zeros(1, J - 1)
    zeros2block = zeros(J - 1, 1)

    hessc = [-2 zeros1block' 2 zeros1block';
             zeros2block -xmat zeros2block xmat;
             2 zeros1block' -2 zeros1block';
             zeros2block xmat zeros2block -xmat]

    hessc = sparse(hessc)

    hess = sparse(f2 + lambda[:eqnonlin][1] * hessc)
    return hess  
end

function mycon(Inputvec, mat, params)
    mid = length(Inputvec) ÷ 2
    Dvec = Inputvec[1:mid]
    Yvec = Inputvec[1+mid:end]
    Pvec = Yvec - Dvec
    ceq = sum(Pvec) - params.Rweight * (Pvec[2:end]' * mat * Pvec[2:end]) - Pvec[1]^2
    return ceq
end


# arbitrarily defined variables
Random.seed!(123)
l_guess = 26
YFmax = 1e-5 .+ (1e-4 - 1e-5) .* rand(13) 
KRshifter = rand(0:0.0001:5.5, 13)

LB = [zeros(13); KRshifter]
UB = [1e3 * ones(13); YFmax .+ KRshifter .+ 1]
guess = [KRshifter; KRshifter .+ 0.001]
struct rParams
    B
end
regionParams = rParams(rand(0:0.009:0.1, 12, 12))
struct Params
    I
    alpha1
    alpha2
    Rweight
end
params = Params(10, 0.3, 0.7, 1.5)
power = rand(0:0.01:0.55, 13, 10)
shifter = rand(0:0.0000001:0.1, 13, 10)
KFshifter = rand(0:0.00001:0.0003, 13)
p_F_LR = 1


# initialize array for x
using JuMP
using Ipopt
x = Vector{Float64}(undef, l_guess)
model = Model(Ipopt.Optimizer)
@variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
@constraint(model, c1, mycon(x, regionParams.B[kk], params) == 0)

# The solve works without the user defined gradient and hessian
@objective(model, Min, obj_f(power, shifter, KFshifter, KRshifter, p_F_LR, params, x))
optimize!(model)



# try to add gradient
x = Vector{AffExpr}(undef, l_guess)
m_g = Model(Ipopt.Optimizer)
@variable(m_g, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])
@constraint(m_g, c1, mycon(x, regionParams.B[kk], params) == 0)

function new_obj_f(x...)
    return obj_f(power, shifter, KFshifter, KRshifter, p_F_LR, params, collect(x))
end
function new_grad_f(g, x...)
    g = grad_f(power, shifter, KFshifter, KRshifter, p_F_LR, params, collect(x))
    return
end
function new_hess_f(h, x...)
    h = f_hess(dual(c1), power, shifter, KFshifter, KRshifter, p_F_LR, regionParams.B, params, x)
    return
end

# if the constraints are user-defined, hessians are disabled
@operator(m_g, op_fg, l_guess, new_obj_f, new_grad_f)
@objective(m_g, Min, op_fg(x...))
optimize!(m_g)


model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])

function new_obj2(x...)
    return obj(power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F_LR, params, collect(x))
end
function new_grad2(g, x...)
    g = grad_f(power[jj], shifter[jj], KFshifter[jj], KRshifter[jj], p_F_LR[jj], params, collect(x))
    return
end
@operator(model, mrkt_eq2, l_guess, new_obj2, new_grad2)
@objective(model, Min, mrkt_eq2(x...))
optimize!(model)


"""
creating an example for interp3 to match up julia and MATLAB
1. replicate MATLABs interp3 using different methods in Julia
the key is to continue using mesh grids over three dimensions

"""

using Interpolations

# Define the coordinate ranges for each dimension
x = 1:0.5:10  # X ranges from 1 to 10 with a step of 0.5
y = 1:0.5:10  # Y ranges from 1 to 10
z = 1:0.5:10  # Z ranges from 1 to 10

# Create a 3D array of data using a simple mathematical function for interpolation
data = [sin(xi) + cos(yi) + zi for xi in x, yi in y, zi in z]

interp = interpolate((collect(x), collect(y), collect(z)), data, Gridded(Linear()))

xi, yi, zi = 5.5, 5.5, 5.5

# Perform the interpolation
interpolated_value = interp(xi, yi, zi)


# ---------------------------------------------------------------------------- #
#                                 tpijs example                                #
# ---------------------------------------------------------------------------- #

using Random
using LinearAlgebra
using Parameters
using Profile


vector = [1.0, 2.0, 3.0]
    # Define a struct to hold parameters
@with_kw struct P
    Q::Int = 3
    tau::Vector{Float64} = vector
    Vs::Matrix{Float64} = rand(3, 4)
    kappa::Float64 = 0.5
    psi::Float64 = 0.7
    Z::Float64 = 1.0
    zsector::Matrix{Float64} = rand(10, Q)
    cdc::Float64 = 1.0
    tau2::Vector{Float64} = tau .+ 1
end


p1 = P()


new = Vector{Float64}(undef, 3)
function test_func!(p1, new::Vector)
    new .= p1.tau .^ p1.psi
    return new
end

@code_warntype test_func!(p1, new)
# Initialize tpijs matrix
tpijs = zeros(Float64, 10, 3)

p = P()

# Compute tpijs for each i
function fill_tpijs!(tpijs::Matrix, p::P, w0::Int, pES::Vector, r::Vector, p_F::Float64)
    for i = 1:p.Q
        @views tpijs[:, i] .= p.tau[i] .* 
                            w0 .^ p.Vs[i, 1] .* 
                            pES .^ (p.Vs[i, 2] .+ p.Vs[i, 3]) .* 
                            (p.kappa .+ (p.kappa .* p_F ./ pES) .^ (1 - p.psi)) .^ 
                            (-(p.psi / (p.psi - 1)) * p.Vs[i, 3]) .*
                            r .^ p.Vs[i, 4] ./
                            (p.Z .* p.zsector[:, i] .* p.cdc)
    end
end

w0 = 2
pES = rand(10)
r = rand(10)
p_F = 0.5
@code_warntype fill_tpijs!(tpijs, p, w0, pES, r, p_F)



#####################
# foo! vs foo() + return

x = Vector{Float64}(undef, 10000)

function foo!(x::Vector)
    x .= rand(10000) .* 0.28
end

function foo()
    t = rand(10000) .* 0.28
    return t
end

BT.@btime foo!(x)
BT.@btime t = foo()



########################################################
using JuMP, Ipopt

# constraint function
function mycon(Inputvec::Vector, B::Matrix, R)
    mid = length(Inputvec) ÷ 2

    Dvec = @view Inputvec[1:mid]
    Yvec = @view Inputvec[1+mid:end]
    Pvec = Yvec .- Dvec

    ceq = sum(Pvec) .- R * (Pvec[2:end]' .* B .* Pvec[2:end]) .- Pvec[1]^2
    return ceq
end

# example data
R = 1.5
l_guess = 10
LB = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0179, 0.0, 0.005, 0.0398]
UB = [1000, 1000, 1000, 1000, 1000, 1000, 1.1223, 0.495, 0.1, 0.126]
guess = [0.0216, 0.0169, -0.001, 0.00377, 0.03879, -0.001, 0.02256, 0.01786, 0.0, 0.004769]
B = 0.005 .+ (0.9 - 0.005) .* rand(4, 4)



x = Vector{Float64}(undef, l_guess)
mid = l_guess ÷ 2
model = Model(Ipopt.Optimizer)
@variable(model, LB[i] <= x[i=1:l_guess] <= UB[i], start=guess[i])

# no longer works in the larger model
# @constraint(model, mycon(x, B, R) .== 0)

# construct the constraint
@expression(model, Pvec, x[1+mid:end] .- x[1:mid])
@expression(model, qmat, -R .* (Pvec[2:end]' .* B .* Pvec[2:end]))
add_to_expression!(qmat, sum(Pvec))


# ---------------------------------------------------------------------------- #
#                              Summary of Model.jl                             #
# ---------------------------------------------------------------------------- #
module Model
export solve_model, data_set_up

using JuMP, Ipopt
import Main.ObjectiveFunction: obj
import Main.Data: data

function data_set_up(data)
    local input = create_model_inputs(data)
    # uses data to create inputs for model
    # inputs include lower bound, upper bound, initial guess, length of initial guess, and 
    # mid (half of the length of the initial guess)

    return input 
end

function add_model_variable(model::Model, input)
    local x = Vector{AffExpr}(undef, input.l_guess)
    for i in eachindex(x)
        local x[i] = AffExpr(0.0)
    end
    return @variable(model, input.LB[i] <= x[i=1:input.l_guess] <= input.UB[i], start = input.guess[i])
end

function add_model_constraint(model::Model, input)
    local x = model[:x]
    local inter = @expression(model, x[input.mid+1:end] .- x[1:input.mid])
    local quadratic = @expression(model, input.quadratic) # results in a 1x1 vector
    add_to_expression!(quadratic[1], inter)

    return @constraint(model, c1, quadratic[1] == 0)
end

function add_model_objective(model::Model, input)
    local x = model[:x]
    return @objective(model, Min, obj(x, input))
end

function solve_model(input)
    local model = Model(Ipopt.Optimizer)
    add_model_variable(model, input)
    add_model_constraint(model, input)
    add_model_objective(model, input)
    optimize!(model)
    return value.(model[:x])

end
end

module Optimize
import Main.Data: data
import Main.Model

# initialize data to update
updated_data = #struct of data to update
chunk = # determine chunk sizes

function optimize_data!(data, updated_data, chunk)
    Threads.@threads for kk in 1:chunk
        local input = data_set_up(data)
        local ouput = solve_model(input)

        updated_data[kk] .= update_data(ouput)
        # the chunks never overlap, so no two threads are simultaneously accessing the same data

    end
end

optimize_data!(data, updated_data, chunk)

module start_w
export w
w = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

end

module shift
import Main.start_w: w
export shift_w, sh_w
function shift_w(w::Vector)
    kk = 1
    while kk < 3
        kk += 1
        w .+= 1
    end

    return w
end
sh_w = shift_w(w)
end



end
