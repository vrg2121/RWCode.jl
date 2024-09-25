"""
To set up the package, in the Julia REPL type the following:

using RWCodeJulia
using RWCodeJulia.ParameterizeRun

# configure model by following prompts
UserConfig = ModelConfig()

# function that runs the model
run_rwcode(UserConfig)

"""


module RWCodeJulia

export ModelConfig, run_rwcode

# ---------------------------------------------------------------------------- #
#                           Prompt User Configuration                          #
# ---------------------------------------------------------------------------- #

include("./ModelConfiguration.jl")
import .ModelConfiguration: ModelConfig
println("Set model configurations by creating struct using ModelConfig(). For example: config = ModelConfig()")
#println("Create path to folder to store results.")

# ---------------------------------------------------------------------------- #
#                         Import all Relevant Functions                        #
# ---------------------------------------------------------------------------- #
# Set Up Parameters
include("./functions/ParamsFunctions.jl") 
include("./Params.jl")
import .Params: setup_parameters

# Load Data
include("./functions/DataAdjustments.jl")
include("./functions/DataLoadsFunc.jl")
include("./DataLoads.jl")
import .DataLoads: load_data

# Initial Equilibrium
include("./functions/MarketEquilibrium.jl")
include("./RegionModel.jl")
include("./functions/MarketFunctions.jl")
include("./Market.jl")
import .Market: solve_market

# Long Run Equilibrium
include("./functions/SteadyStateFunctions.jl")
include("./SteadyState.jl")
import .SteadyState: solve_steadystate

# Run Transitional Dynamics
include("./functions/TransitionFunctions.jl")
include("./Transition.jl")
import .Transition: solve_transition

# Write Data
include("./WriteData.jl")
import .WriteData: writedata

# Write Data with Battery configurations
include("./WriteDataBattery.jl")
import .WriteDataBattery: writedata_battery

# -------------------- Run Exogenous Technology Equilibria ------------------- #

# Long Run Equilibrium with Exogenous Tech
include("./functions/SteadyStateExogFunc.jl")
include("./SteadyStateExog.jl")
import .SteadyStateExog: solve_steadystate_exog

# Transitional Dynamics with Exogenous Tech
include("./functions/TransitionExogFunc.jl")
include("./TransitionExog.jl")
import .TransitionExog: solve_transition_exog

# Data Outputs with Exogenous Tech
include("./WriteDataExog.jl")
import .WriteDataExog: writedata_exog


function run_rwcode(config::ModelConfig)

    # ---------------------------------------------------------------------------- #
    #                               Set Up Parameters                              #
    # ---------------------------------------------------------------------------- #

    println("Setting up parameters...")
    P = setup_parameters();

    # ---------------------------------------------------------------------------- #
    #                                  Data Loads                                  #
    # ---------------------------------------------------------------------------- #
    
    println("Loading data inputs...")
    D = load_data(P)
    
    # ---------------------------------------------------------------------------- #
    #                              Initial Equilibrium                             #
    # ---------------------------------------------------------------------------- #

    println("Solving initial equilibrium...")    
    M = solve_market(P, D, config);
    
    
    # ---------------------------------------------------------------------------- #
    #                             Long Run Equilibrium                             #
    # ---------------------------------------------------------------------------- #

    println("Solving initial long run equilibrium...")
    S = solve_steadystate(P, D, M, config)

    # ---------------------------------------------------------------------------- #
    #                           Run Transitional Dynamics                          #
    # ---------------------------------------------------------------------------- #
    if config.RunTransition == 1

        # ---------------------- Run Transition without Subsidy ---------------------- #
        println("Solving transitional dynamics without Subsidy...")
        Subsidy = 0
        T = solve_transition(P, D, M, S, Subsidy, config)

        println("Writing outputs of transitional dynamics without Subsidy...")
        writedata(P, D, M, S, T, Subsidy, config)
        #writedata(P, D, M, S, T, Subsidy, config, path)

        # ------------------------ Run Transition with Subsidy ----------------------- #
        println("Solving transitional dynamics with Subsidy...")
        Subsidy = 1
        TS = solve_transition(P, D, M, S, Subsidy, config)

        println("Writing outputs of transitional dynamics with Subsidy...")
        writedata(P, D, M, S, TS, Subsidy, config)
        #writedata(P, D, M, S, TS, Subsidy, config, path)

    end

    # ---------------------------------------------------------------------------- #
    #                              Battery Robustness                              #
    # ---------------------------------------------------------------------------- #
    """
    If you run config.RunTransition == 1, data will be written using writedata without Subsidy and then with Subsidy.
    If you run config.RunTransition == 1 && config.RunBatteries == 1, data will be written 
        using writedata (with subsidy, without subsidy) and then written using writedatabattery.
    When you run config.RunTransition == 1, should hoursofstorage always be 0?
    """
    
    if config.RunBatteries == 1

        for bb in 1:length(config.hoursvec)            
            config.hoursofstorage = config.hoursvec[bb]
            Subsidy = 0

            println("Solving long run equilibrium when battery storage hours = $(config.hoursofstorage)...")
            SB = solve_steadystate(P, D, M, config)

            println("Solving transitional dynamics when battery storage hours = $(config.hoursofstorage)...")
            TB = solve_transition(P, D, M, SB, Subsidy, config)
            
            println("Writing outputs when battery storage hours = $(config.hoursofstorage)...")
            writedata_battery(P, M, SB, TB, config)
        end
    end

    # ---------------------------------------------------------------------------- #
    #                               Grid Improvements                              #
    # ---------------------------------------------------------------------------- #
    if config.RunImprovement == 1
        #include("./SteadyStateImp.jl")
        #SI = solve_steadystate_imp()
        #writedata_imp(P, D, M, SI, config)
    end

    # ---------------------------------------------------------------------------- #
    #                                Exogenous Tech                                #
    # ---------------------------------------------------------------------------- #
    if config.RunExog == 1
        config.hoursofstorage = 0
        for exogindex in 3:-1:1
            println("Solving long run equilibrium with exogenous tech when exog index = $exogindex...")
            SE = solve_steadystate_exog(P, D, M, config, exogindex)

            println("Solving transitional dynamics with exogenous tech when exog index = $exogindex...")
            TE = solve_transition_exog(P, D, M, SE, config, exogindex)
            
            println("Writing outputs with exogenous tech when exog index = $exogindex...")
            writedata_exog(TE, exogindex)
        end
            
    end

    # ---------------------------------------------------------------------------- #
    #                                  Curtailment                                 #
    # ---------------------------------------------------------------------------- #
    """
    User defined hours of storage is overwritten when RunCurtailment == 1.
    """
    if config.RunCurtailment == 1
        Subsidy = 0
        config.hoursofstorage = 12

        println("Solving long run equilibrium with curtailment (hours of battery storage = 12)...")
        SC = solve_steadystate(P, D, M, config)

        println("Solving transitional dynamics with curtailment (hours of battery storage = 12)...")
        TC = solve_transition(P, D, M, SC, Subsidy, config)

        println("Writing output for curtailment (hours of battery storage = 12)...")
        writedata_battery(P, M, SC, TC, config)
    end


end

end
