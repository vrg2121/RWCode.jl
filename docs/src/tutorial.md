# Tutorial: Run the Model

Start this tutorial after you have completed the [Installation](@ref).

#### Table of Contents
```@contents
Pages = ["tutorial.md"]
Depth = 2:3
```

## PC and Git
Before running the model, check that the following are complete from the installation:

1. Data and Guesses and the RWCode Repository have been downloaded.

1. The Results folder has been created.

1. Julia version is 1.10

1. There are 12+ threads in the Julia session.

Activate and compile the package:
- Navigate to the working directory with the package...

```julia
julia> cd("PATH/TO/RWCode")

julia> ]

(@v1.10) pkg> activate .

(RWCode) pkg> instantiate

julia> using RWCode
```

Set paths for Data, Guesses and Results Folders:

```julia
julia> D = "FULL/PATH/TO/DATA"
"FULL/PATH/TO/DATA"

julia> G = "FULL/PATH/TO/Guesses"
"FULL/PATH/TO/Guesses"

julia> R = "FULL/PATH/TO/Results"
"FULL/PATH/TO/Results"
```

Configure the model (see Model Configuration for more details):

```julia
julia> config = ModelConfig()
```

Run the model:

```julia
julia> run_rwcode(config, D, G, R)
```


## Docker
1. There are 12+ threads in the Julia session.

Set up variable paths to data, guesses and results:

```julia
julia> D = "/home/jl/Data"

julia> G = "/home/jl/Guesses"

julia> R = "/home/jl/Results"
```

!!! note 
	The location of the Data, Guesses and Results files are pre-set when using the Docker image.

Configure the model:

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
```
* In this example, all of the values are automatically set to the default. For more information about the model configurations see this page.


Run the model:

```julia
julia> run_rwcode(config, D, G, R)
```

!!! info "Locating Outputs"
	All of the outputs will be labeled and stored as .csv files in the Results file at `/home/jl/Results`.

Copy the results from container to local:
* In the terminal/CLI...

```bash
docker ps -a	# this command will print the container ID in the first column.
docker cp <containerID>:/home/jl/Results </host/path/target>
```

* Change `</host/path/target>` to the local path you want to copy the data into.
