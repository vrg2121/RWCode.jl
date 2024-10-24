# Model Configuration
This model has various configurations that the user can set using the function:

```julia
julia> config = ModelConfig()
```

## Input
This function guides the user through creating a mutable struct with the model configurations. The outputs for this function should be stored in a variable, e.g. `config = ModelConfig()`. 

  * `RunTransition`: Binary Input, i.e.`(1, 0)`

      + When set to 1, the model will run a 20 year transition from 2020-2040 as wind and solar are incorporated into the global energy grid.
      + The model automatically runs the transition with a subsidy included in the prices for the US.

  * `RunBatteries`: Binary Input, i.e. `(1, 0)`

      + When set to 1, the model will run the 20 year transition with the hours of storage for batteries set by the user in hoursvec.
      + `RunBatteries = 1` and `hoursvec` must both be set.
      + By default `hoursvec = [2, 4, 6]`

  * `RunExog`: Binary Input, i.e. `(1, 0)`

      + When set to 1, the model assumes that technology is exogenous. 
      + Calculations for the steady state equilibrium and the transition assume that technology does not change.
      + Specifically, the exogenous values impact the guesses for renewable capital prices in both the long run and the transition, setting the guesses on projections for wind / solar respectively.

  * `RunCurtailment`: Binary Input, i.e. `(1, 0)`

      + When set to 1, there is no subsidy during transition and battery hours of storage is assumed to be 12.
      + This problem demonstrates the effects of renewable intermittency (inconsitency in energy production) without optimized storage hours.

  * `Transiter`: Integer Input, i.e. `100`

      + Controls the number of iterations the transition model will run.
      + To create data for the paper, the transition model runs for 100 iterations, i.e. `Transiter = 100`.
      + To check that the code is working Transiter = 2.

  * `Initialprod`: Integer Input, i.e. `100`

      + Initial level of renewables cost of production.
      + To create data for the paper, `Initialprod = 100`.

  * `hoursofstorage`: Integer Input, i.e. `4`

      + Hours of battery storage for renewable energy.
      + This value is always set by the model. User input will be overwritten.
      + Use `hoursvec` to set battery hours of prices when `RunBatteries==1`.

  * `hoursvec`: List Input, i.e. `2,4,6`

      + The hoursvec is a vector of battery storage hours used as inputs for the model. 
      + This will output results for the change in battery prices over the transition given that the number of hours of battery storage. 
      + It is useful for comparing different hours of storage.
      + This variable is only relevant when `RunBatteries==1`.


## Examples
### Default Configuration
```julia
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
```


### Run Batteries with Hours 1,3,5
```julia
julia> config = ModelConfig()
Enter RunTransition (0 or 1, default = 1):
0
Enter RunBatteries (0 or 1, default=0):
1
Enter RunExog (0 or 1, default=0):
0
Enter RunCurtailment (0 or 1, default=0):
0
Enter the Number of Transition Iterations (recommend 0-100 iterations, default=2):
2
Enter Initial Production (default = 100):
0
Enter hoursofstorage (default=0):
0
Enter hoursvec (comma-separated, default = 2,4,6):
1,3,5
ModelConfig(0, 1, 0, 0, 2, 0, 0, [1, 3, 5])
```