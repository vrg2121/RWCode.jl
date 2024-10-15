"""
To set up the package:

# ---------------------------------------------------------------------------- #
#        Clone the Package and Set Environment with Package Dependencies       #
# ---------------------------------------------------------------------------- #

# clone the package from git
shell> git clone https://github.com/vrg2121/RWCode.jl.git
Cloning into 'RWCode.jl'...
...

# go to your preferred julia IDE or in the CLI simply type "julia"
julia> ]

(@v1.10) pkg> activate RWCode.jl
Activating project at `~/RWCode.jl`

(RWCode) pkg> instantiate

# for a detailed explanation see: https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project


# ---------------------------------------------------------------------------- #
#                    Download Data and Guesses for the Model                   #
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
#              Set Paths to Data and Guesses, Configure Model Run              #
# ---------------------------------------------------------------------------- #

using RWCodeJulia
using RWCodeJulia.ModelConfiguration

# configure model by following prompts
config = ModelConfig()

# set path to the Data folder
D = "C:/Users/UserName/path/to/Data"

# set path to the Guesses folder
G = "C:/Users/UserName/path/to/Guesses"

# set path for Results output
R = "C:/Users/UserName/path/to/Results"

# function that runs the model
run_rwcode(config, D, G, R)

"""


module RWCode

export ModelConfig, run_rwcode

# ---------------------------------------------------------------------------- #
#                           Prompt User Configuration                          #
# ---------------------------------------------------------------------------- #

include("./ModelConfiguration.jl")
import .ModelConfiguration: ModelConfig
# Set model configurations by creating struct using ModelConfig(). For example:
## config = ModelConfig()

######################################################
# create path to Results folder
# R = "C:/Users/vrg2121/.julia/dev/RWCodeJulia/Results"

######################################################
# create path to Data folder

# D = "C:/Users/vrg2121/.julia/dev/RWCodeJulia/Data"

#######################################################
# create path to Guesses folder
# G = "C:/Users/vrg2121/.julia/dev/RWCodeJulia/Guesses"



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

"""
    run_rwcode(config::ModelConfig, Data::String, Guesses::String, Results::String)

Takes in user configuration, data, and guesses to run the model. Output is saved in Results file.

## Input
- `config::ModelConfig`-- struct of user defined model configurations.
- `Data::String`-- full path to Data folder.
- `Guesses::String`-- full path to Guesses folder.
- `Results::String`-- full path to Results folder.

## Output
CSV files in Results folder. Name and content of files depend on user configuration.
"""
function run_rwcode(config::ModelConfig, D::String, G::String, R::String)

    # ---------------------------------------------------------------------------- #
    #                               Set Up Parameters                              #
    # ---------------------------------------------------------------------------- #

    println("Setting up parameters...")
    P = setup_parameters(D, G);

    # ---------------------------------------------------------------------------- #
    #                                  Data Loads                                  #
    # ---------------------------------------------------------------------------- #
    
    println("Loading data inputs...")
    DL = load_data(P, D)
    
    # ---------------------------------------------------------------------------- #
    #                              Initial Equilibrium                             #
    # ---------------------------------------------------------------------------- #

    println("Solving initial equilibrium...")    
    M = solve_market(P, DL, config, G);
    
    
    # ---------------------------------------------------------------------------- #
    #                             Long Run Equilibrium                             #
    # ---------------------------------------------------------------------------- #

    println("Solving initial long run equilibrium...")
    S = solve_steadystate(P, DL, M, config, G)

    # ---------------------------------------------------------------------------- #
    #                           Run Transitional Dynamics                          #
    # ---------------------------------------------------------------------------- #
    if config.RunTransition == 1

        # ---------------------- Run Transition without Subsidy ---------------------- #
        println("Solving transitional dynamics without Subsidy...")
        Subsidy = 0
        T = solve_transition(P, DL, M, S, Subsidy, config, G)

        println("Writing outputs of transitional dynamics without Subsidy...")
        writedata(P, DL, M, S, T, Subsidy, config, R)
        #writedata(P, DL, M, S, T, Subsidy, config, path)

        # ------------------------ Run Transition with Subsidy ----------------------- #
        println("Solving transitional dynamics with Subsidy...")
        Subsidy = 1
        TS = solve_transition(P, DL, M, S, Subsidy, config, G)

        println("Writing outputs of transitional dynamics with Subsidy...")
        writedata(P, DL, M, S, TS, Subsidy, config, R)
        #writedata(P, DL, M, S, TS, Subsidy, config, path)

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
            SB = solve_steadystate(P, DL, M, config, G)

            println("Solving transitional dynamics when battery storage hours = $(config.hoursofstorage)...")
            TB = solve_transition(P, DL, M, SB, Subsidy, config, G)
            
            println("Writing outputs when battery storage hours = $(config.hoursofstorage)...")
            writedata_battery(P, M, SB, TB, config, R)
        end
    end


    # ---------------------------------------------------------------------------- #
    #                                Exogenous Tech                                #
    # ---------------------------------------------------------------------------- #
    if config.RunExog == 1
        config.hoursofstorage = 0
        for exogindex in 3:-1:2
            # exogindex == 1 does not converge
            println("Solving long run equilibrium with exogenous tech when exog index = $exogindex...")
            SE = SteadyStateExog.solve_steadystate_exog(P, DL, M, config, exogindex, G)


            println("Solving transitional dynamics with exogenous tech when exog index = $exogindex...")
            TE = solve_transition_exog(P, DL, M, SE, config, exogindex)
            
            println("Writing outputs with exogenous tech when exog index = $exogindex...")
            writedata_exog(TE, exogindex, R)
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
        SC = solve_steadystate(P, DL, M, config, G)

        println("Solving transitional dynamics with curtailment (hours of battery storage = 12)...")
        TC = solve_transition(P, DL, M, SC, Subsidy, config, G)

        println("Writing output for curtailment (hours of battery storage = 12)...")
        writedata_battery(P, M, SC, TC, config, R)
    end


end

end
