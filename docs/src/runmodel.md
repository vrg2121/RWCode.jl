# Run Model with run_rwcode()
Users can run the model using the function:
```julia
julia> run_model(config, D, G, R)
```

This function runs the model with the user Model Configuration and places outputs in the `Results` folder.

```julia
run_rwcode(config::ModelConfig, Data::String, Guesses::String, Results::String)
```

Takes in user configuration, data, and guesses to run the model. Output is saved in Results file.

## Input
  * `config::ModelConfig`-- struct of user defined model configurations.
  * `Data::String`-- full path to Data folder.
  * `Guesses::String`-- full path to Guesses folder.
  * `Results::String`-- full path to Results folder.

## Output
CSV files in `Results` folder. Name and content of files depend on user configuration.

## Example
```julia
julia> D = "PATH/TO/DATA"
"PATH/TO/DATA"

julia> G = "PATH/TO/GUESSES"
"PATH/TO/GUESSES"

julia> R = "PATH/TO/RESULTS"
"PATH/TO/RESULTS"

julia> config = ModelConfig()
Enter RunTransition (0 or 1, default = 1):

Enter RunBatteries (0 or 1, default=0):

Enter RunExog (0 or 1, default=0):

Enter RunCurtailment (0 or 1, default=0):

Enter the Number of Transition Iterations (recommend 0-100 iterations, default=2):

Enter Initial Production (default = 100):

Enter hoursofstorage (default=0):

Enter hoursvec (comma-separated, default = 2,4,6):

ModelConfig(1, 0, 0, 0, 2, 100, 0, [2, 4, 6])

julia> run_model(config, D, G, R)
```

