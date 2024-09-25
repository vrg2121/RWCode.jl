# Parameterize Run
module ModelConfiguration

export ModelConfig

function prompt_for_bin(prompt, default=nothing, binary=false)
    while true
        println(prompt)
        val = readline()
        if isempty(val) && default !== nothing
            return default
        else
            try
                parsed_val = parse(Int, val)
                if binary && !(parsed_val in (0, 1))
                    println("Error: Only 0 or 1 is allowed. Please try again.")
                else
                    return parsed_val
                end
            catch
                println("Invalid input. Please enter a valid integer.")
            end
        end
    end
end

function prompt_for_int(prompt, default=nothing)
    println(prompt)
    val = readline()
    return isempty(val) && default !== nothing ? default : parse(Int, val)
end

function prompt_for_vector(prompt, default=nothing)
    println(prompt)
    val = readline()
    if isempty(val) && default !== nothing
        return default
    else
        #val = map(x -> parse(Int64, x), val)
        return parse.(Float64, split(val, ","))
    end
end


mutable struct ModelConfig
    RunTransition::Int64
    RunBatteries::Int64
    RunImprovement::Int64 
    RunExog::Int64
    RunCurtailment::Int64
    Transiter::Int64
    Initialprod::Int64
    hoursofstorage::Int64
    hoursvec::Vector{Int64}

    function ModelConfig()
        RunTransition = prompt_for_bin("Enter RunTransition (0 or 1, default = 1):", 1, true)
        RunBatteries = prompt_for_bin("Enter RunBatteries (0 or 1, default=0):", 0, true)
        RunImprovement = prompt_for_bin("Enter RunImprovement (0 or 1, default=0):",0, true)
        RunExog = prompt_for_bin("Enter RunExog (0 or 1, default=0):",0, true)
        RunCurtailment = prompt_for_bin("Enter RunCurtailment (0 or 1, default=0):",0, true)
        Transiter = prompt_for_int("Enter the Number of Transition Iterations (recommend 0-100 iterations, default=2):", 2)
        Initialprod = prompt_for_int("Enter Initial Production (default = 100):", 100)
        hoursofstorage = prompt_for_int("Enter hoursofstorage (default=0):", 0)
        hoursvec = prompt_for_vector("Enter hoursvec (comma-separated, default = 2,4,6):", [2.0, 4.0, 6.0])
        
        new(RunTransition, RunBatteries, RunImprovement, RunExog, RunCurtailment, Transiter, Initialprod, hoursofstorage, hoursvec)
    end
    
end


end